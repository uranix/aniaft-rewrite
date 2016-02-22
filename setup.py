from distutils.core import setup, Extension

setup(
    name = 'aniaft',
    version = '0.1',
    ext_modules = [
        Extension("_aniaft",
            sources=["boundary.cpp", "tree.cpp", "tria.cpp", "mesh.cpp", "aniaft.i"],
            extra_compile_args=['-std=c++11'],
            swig_opts=['-modern', '-c++']
        )
    ]
)
