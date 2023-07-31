fn main() {
    println!("cargo:rustc-link-lib=dylib=Ole32");
    println!("cargo:rustc-link-lib=dylib=Shell32");
}
