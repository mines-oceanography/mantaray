"""
How to use it

Compile the Rust library into ./tmp
cargo cinstall --manifest-path=Cargo.toml --destdir=./tmp --prefix=ray --libdir=ray

Run Python that will import the library and use it
python demo.py
"""
import ctypes

# dylib or so?
rust = ctypes.CDLL("tmp/ray/libray_tracing.dylib")

rust.single_ray.argtypes = (
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
)
# rust.single_ray.restype = ctypes.c_double


if __name__ == "__main__":
    rust.single_ray(
        "./island.nc".encode("utf-8"),
        -1000.0, 0.0, 0.01, 0.0, 10, 2)
