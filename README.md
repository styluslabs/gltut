# gltut #

This is a project that was used for learning OpenGL and exploring techniques for rendering vector graphics in OpenGL.

`lbmain.cpp` is the `main()` file (gltut.cpp and main.cpp are older explorations) - see `main()` for the accepted keyboard input.

It uses [stb_truetype](https://github.com/nothings/stb) to load a font and then renders text with one of several different techniques, implemented in both OpenGL (`gl_*.cpp`) and with a software renderer (`sw_*.cpp`).  The software renderer tries to implement a kind of GPU emulator (`swrender.cpp`), such that GLSL shaders can be used with minimial modification as C code.  It does not currently implement anything related to depth, stencil, etc. (but it would be interesting to add this).

Included are signed distance field (SDF) and various "exact coverage" techniques, some of which were ultimately used for [nanovgXC](https://github.com/styluslabs/nanovgXC).


### Dependencies ###

* [stb](https://github.com/nothings/stb) libraries and [glm](https://glm.g-truc.net/) are assumed to be in the parent of this directory; edit `Makefile` as needed if these libraries are in a different location
* requires [SDL2](https://www.libsdl.org/) and [GLEW](https://github.com/nigels-com/glew/releases).  On Debian/Ubuntu, just install the `libsdl2-dev` and `libglew-dev` packages.  On Windows, open Makefile.msvc to see or change the expected locations for the libraries and headers.


### Building ###

After installing dependencies:

Linux:
* run `make`

Windows:
* requires Visual C++ (free Visual Studio Community Edition is fine)
* download GNU make built for Windows: http://www.equation.com/servlet/equation.cmd?fa=make and ensure it is available in the path
* open a Visual Studio command prompt (from the Start Menu) and run `make`


### Credits ###

* https://github.com/hansent/lbfont
