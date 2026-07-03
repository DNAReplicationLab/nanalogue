const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize: std.builtin.OptimizeMode = .ReleaseSafe;

    const exe = b.addExecutable(.{
        .name = "reads_table",
        .root_module = b.createModule(.{
            .root_source_file = b.path("reads_table.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    exe.root_module.addIncludePath(.{ .cwd_relative = "htslib_compiled/include" });
    exe.root_module.addLibraryPath(.{ .cwd_relative = "htslib_compiled/lib" });

    exe.root_module.linkSystemLibrary("hts", .{});
    exe.root_module.linkSystemLibrary("z", .{});
    exe.root_module.linkSystemLibrary("bz2", .{});
    exe.root_module.linkSystemLibrary("lzma", .{});
    exe.root_module.linkSystemLibrary("deflate", .{});
    exe.root_module.linkSystemLibrary("curl", .{});
    exe.root_module.linkSystemLibrary("crypto", .{});
    exe.root_module.linkSystemLibrary("pthread", .{});

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| run_cmd.addArgs(args);

    const run_step = b.step("run", "Run reads_table");
    run_step.dependOn(&run_cmd.step);
}
