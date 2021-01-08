from distutils.core import setup, Extension

# http://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module


from distutils.command.build import build

class CustomBuild(build):
    sub_commands = [
        ('build_ext', build.has_ext_modules),
        ('build_py', build.has_pure_modules),
        ('build_clib', build.has_c_libraries),
        ('build_scripts', build.has_scripts),
    ]

setup(
    name = 'aniaft',
    version = '0.2',
    cmdclass = { 'build' : CustomBuild },
    ext_modules = [
        Extension("_aniaft",
            sources=["boundary.cpp", "tree.cpp", "tria.cpp", "mesh.cpp", "aniaft.i"],
            extra_compile_args=['-std=c++11', '-Wall'],
            swig_opts=['-c++', '-Wall']
        )
    ],
    py_modules = ['aniaft']
)
