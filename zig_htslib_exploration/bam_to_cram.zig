const std = @import("std");
const c = @cImport({
    @cInclude("htslib/sam.h");
    @cInclude("htslib/hts.h");
    @cInclude("sys/stat.h");
});

// Prefix view needed because C bitfields make bam1_t opaque to Zig.
const Bam1CoreView = extern struct {
    core: c.bam1_core_t,
    id: u64,
    data: [*]u8,
    l_data: i32,
    m_data: u32,
    _unused_bitfield_storage: u32,
};

// Check if paths are same by retrieving file system metadata if available.
fn isPathsSame(first: [*:0]const u8, second: [*:0]const u8) bool {
    var first_stat: c.struct_stat = undefined;
    var second_stat: c.struct_stat = undefined;
    if (c.stat(first, &first_stat) != 0 or c.stat(second, &second_stat) != 0)
        return false;
    return first_stat.st_dev == second_stat.st_dev and first_stat.st_ino == second_stat.st_ino;
}

pub fn main(init: std.process.Init) !void {
    // Allocate memory
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);

    // Process command line arguments and check for duplicate file paths.
    if (args.len != 4 and args.len != 5) {
        std.debug.print(
            "usage: {s} <input.bam> <output.cram> <reference.fa> [threads]\n",
            .{args[0]},
        );
        std.process.exit(1);
    }
    const index_path = try std.fmt.allocPrintSentinel(arena, "{s}.crai", .{args[2]}, 0);
    if (std.mem.eql(u8, args[1], args[2]) or
        std.mem.eql(u8, args[1], args[3]) or
        std.mem.eql(u8, args[2], args[3]) or
        std.mem.eql(u8, args[1], index_path) or
        std.mem.eql(u8, args[3], index_path) or
        isPathsSame(args[1].ptr, args[2].ptr) or
        isPathsSame(args[1].ptr, args[3].ptr) or
        isPathsSame(args[2].ptr, args[3].ptr) or
        isPathsSame(args[1].ptr, index_path.ptr) or
        isPathsSame(args[3].ptr, index_path.ptr) or
        isPathsSame(args[2].ptr, index_path.ptr))
        return error.OutputAliasesInput;
    const threads = if (args.len == 5)
        try std.fmt.parseInt(c_int, args[4], 10)
    else
        1;
    if (threads < 1) return error.InvalidThreadCount;

    // Open input
    const input = c.sam_open(args[1].ptr, "r") orelse return error.InputOpenFailed;
    defer _ = c.sam_close(input);

    // Read header
    const header = c.sam_hdr_read(input) orelse return error.HeaderReadFailed;
    defer c.sam_hdr_destroy(header);

    // Open output
    const output = c.sam_open(args[2].ptr, "wc") orelse return error.OutputOpenFailed;
    var output_is_open = true;
    defer if (output_is_open) {
        _ = c.sam_close(output);
    };

    // Set options, write header, and start index
    if (c.hts_set_opt(output, c.CRAM_OPT_VERSION, "3.1") != 0)
        return error.CramVersionFailed;
    if (c.hts_set_opt(output, c.CRAM_OPT_EMBED_REF, @as(c_int, 0)) != 0 or
        c.hts_set_opt(output, c.CRAM_OPT_NO_REF, @as(c_int, 0)) != 0)
        return error.ExternalReferenceSetupFailed;
    if (c.hts_set_opt(output, c.CRAM_OPT_STORE_MD, @as(c_int, 1)) != 0 or
        c.hts_set_opt(output, c.CRAM_OPT_STORE_NM, @as(c_int, 1)) != 0)
        return error.TagPreservationSetupFailed;
    if (c.hts_set_fai_filename(output, args[3].ptr) != 0)
        return error.ReferenceLoadFailed;
    if (c.hts_set_threads(output, threads) != 0)
        return error.ThreadSetupFailed;
    if (c.sam_hdr_write(output, header) != 0)
        return error.HeaderWriteFailed;
    if (c.sam_idx_init(output, header, 0, index_path.ptr) != 0)
        return error.IndexInitFailed;

    // Initialize record
    const record = c.bam_init1() orelse return error.RecordAllocFailed;
    defer c.bam_destroy1(record);

    // Loop through records, writing CRAM
    var count: u32 = 0;
    var status: c_int = -1; // -1 means end of stream
    var previous_tid: i32 = -1;
    var previous_pos: i64 = -1;
    var saw_unplaced = false;
    while (true) {
        // Impose a ceiling on number of records
        if (count < std.math.maxInt(u32)) {
            count += 1;
        } else {
            return error.TooManyRecords;
        }

        // Read data
        status = c.sam_read1(input, header, record);
        if (status < 0) break;

        // Check if data is in order
        const record_view: *align(1) const Bam1CoreView = @ptrCast(record);
        const tid = record_view.core.tid;
        const pos = record_view.core.pos;
        if (tid < 0) {
            saw_unplaced = true;
            if (pos > -1) return error.InvalidUnmappedRead;
        } else {
            if (pos < 0) return error.InvalidMappedRead;
            if (saw_unplaced or ((tid < previous_tid or (tid == previous_tid and pos < previous_pos))))
                return error.RecordsNotCoordinateSorted;
            previous_tid = tid;
            previous_pos = pos;
        }

        // Error if writing data fails
        if (c.sam_write1(output, header, record) < 0)
            return error.RecordWriteFailed;
    }
    // Return if errors in reading data or writing index data
    if (status < -1) return error.RecordReadFailed;
    if (c.sam_idx_save(output) < 0) return error.IndexSaveFailed;

    // Close output
    const close_status = c.sam_close(output);
    output_is_open = false;
    if (close_status < 0) return error.OutputCloseFailed;
}
