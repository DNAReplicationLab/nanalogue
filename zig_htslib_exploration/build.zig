const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize: std.builtin.OptimizeMode = .ReleaseSafe;

    const reads_table = b.addExecutable(.{
        .name = "reads_table",
        .root_module = b.createModule(.{
            .root_source_file = b.path("reads_table.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    const bam_to_cram = b.addExecutable(.{
        .name = "bam_to_cram",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bam_to_cram.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    const executables = [_]*std.Build.Step.Compile{ reads_table, bam_to_cram };
    for (executables) |exe| {
        configureHtslib(exe);
        b.installArtifact(exe);
    }
}

fn configureHtslib(exe: *std.Build.Step.Compile) void {
    exe.root_module.addIncludePath(.{ .cwd_relative = "htslib_compiled/include" });

    exe.root_module.addObjectFile(.{ .cwd_relative = "htslib_compiled/lib/libhts.a" });
    exe.root_module.linkSystemLibrary("z", .{});
    exe.root_module.linkSystemLibrary("bz2", .{});
    exe.root_module.linkSystemLibrary("lzma", .{});
    exe.root_module.linkSystemLibrary("deflate", .{});
    exe.root_module.linkSystemLibrary("curl", .{});
    exe.root_module.linkSystemLibrary("crypto", .{});
    exe.root_module.linkSystemLibrary("pthread", .{});
}
