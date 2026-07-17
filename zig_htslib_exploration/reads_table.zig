const std = @import("std");
const c = @cImport({
    @cInclude("htslib/sam.h");
    @cInclude("htslib/hts.h");
});
const assert = std.debug.assert;

// Prefix view of htslib's bam1_t.
// We use this because @cImport makes bam1_t opaque due to bitfields.
// This must match the leading layout of bam1_t in the htslib version we build against.
// TODO: Add a comptime check that bam1_t's definition did not change in htslib.
const Bam1CoreView = extern struct {
    core: c.bam1_core_t,
    id: u64,
    data: [*]u8,
    l_data: i32,
    m_data: u32,
    _unused_bitfield_storage: u32,
};

const Alignment_types = enum {
    primary_forward,
    primary_reverse,
    secondary_forward,
    secondary_reverse,
    supplementary_forward,
    supplementary_reverse,
    unmapped,

    pub fn name(self: Alignment_types) []const u8 {
        return switch (self) {
            Alignment_types.primary_forward => "primary_forward",
            Alignment_types.primary_reverse => "primary_reverse",
            Alignment_types.secondary_forward => "secondary_forward",
            Alignment_types.secondary_reverse => "secondary_reverse",
            Alignment_types.supplementary_forward => "supplementary_forward",
            Alignment_types.supplementary_reverse => "supplementary_reverse",
            Alignment_types.unmapped => "unmapped",
        };
    }
};

pub fn set_alignment(flag: u16) !Alignment_types {
    return switch (flag) {
        0 => Alignment_types.primary_forward,
        4 => Alignment_types.unmapped,
        16 => Alignment_types.primary_reverse,
        256 => Alignment_types.secondary_forward,
        272 => Alignment_types.secondary_reverse,
        2048 => Alignment_types.supplementary_forward,
        2064 => Alignment_types.supplementary_reverse,
        else => error.AlignmentTypeError,
    };
}

pub fn main(init: std.process.Init) !void {

    // do not want htslib logging to print messages
    c.hts_set_log_level(c.HTS_LOG_OFF);

    var stdout_buffer: [64 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writer(init.io, &stdout_buffer);
    const stdout = &stdout_writer.interface;

    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);

    if (args.len != 2) {
        std.debug.print("usage: {s} <file.bam>\n", .{args[0]});
        std.process.exit(1);
    }

    const path = args[1];

    const fp = c.sam_open(path.ptr, "r");
    if (fp == null) return error.OpenFailed;
    defer _ = c.sam_close(fp);

    const hdr = c.sam_hdr_read(fp);
    if (hdr == null) return error.HeaderReadFailed;
    defer c.sam_hdr_destroy(hdr);

    const rec = c.bam_init1() orelse return error.RecordAllocFailed;
    defer c.bam_destroy1(rec);

    var count: u64 = 0;
    var align_len: i64 = 0;
    var seq_len: i32 = 0;

    try stdout.print("{s}\n", .{"read_id\talign_length\tsequence_length_template\talignment_type"});

    while (c.sam_read1(fp, hdr, rec) >= 0) {
        const rec_view: *align(1) const Bam1CoreView = @ptrCast(rec);
        align_len = if ((rec_view.core.flag & c.BAM_FUNMAP) != 0)
            0
        else
            c.bam_endpos(rec) - rec_view.core.pos;
        seq_len = rec_view.core.l_qseq;

        // simple assertions for validity
        assert(align_len >= 0);
        assert(seq_len >= 0);
        assert(rec_view.core.l_qname > rec_view.core.l_extranul + 1);

        const qname_len = rec_view.core.l_qname - rec_view.core.l_extranul - 1;
        const read_id = rec_view.data[0..qname_len];

        // increment counter and check it is not too large
        count += 1;
        if (count > std.math.maxInt(u32)) return error.TooManyRecords;

        try stdout.print("{s}\t{d}\t{d}\t{s}\n", .{ read_id, align_len, seq_len, blk: {
            const b = set_alignment(rec_view.core.flag) catch |err| {
                return err;
            };
            break :blk b.name();
        } });
    }

    try stdout.flush();
    assert(count > 0);
}
