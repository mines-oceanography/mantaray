#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

int32_t single_ray(const char *bathymetry_path,
                   double x0,
                   double y0,
                   double kx0,
                   double ky0,
                   double end_time,
                   double step_size);

#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus
