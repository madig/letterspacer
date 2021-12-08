profile:
    cargo build --release
    perf record --call-graph dwarf -F 999 target/release/letterspacer ../../DejaVu-Font-Sandbox/sources/cubic/DejaVuSans.ufo
    perf script -F +pid > /tmp/test.perf
