# 3D Ray Tracer

## Compiling and Running

A `Makefile` is provided, so code can be compiled and run with

```bash
make
./raytracer.exe [input_file]
```

Additionally, a `CMakeLists.txt` file is also provided so `CMake` can be used as well.

## Input

The input file **must** contain a line of the form

```text
imsize [width] [height]
```

where `[width]` and `[height]` are both positive integers indicated the desired width and height of
the output image.

There are several other keywords which are supported in the input. Look at
`SampleImage.input` for an example of the input file format.
