
# Todo

- [ ] find meshlink_platform and replace lib search in configure.ac
- [ ] initial mesh via API, not xml parse

# Mac

git clone git@github.com:pointwise/MeshLink.git

cd MeshLink/external
tar xzvf [pw-ge]/ge_lite.macx-clang.6.8.tgz
ln -s ge_lite.macx-clang.6.8/ gelite
cd gelite
ln -s macx-clang macosx

cd MeshLink/src/mlparser_xerces
make -f Makefile.macosx BUILD=Release

cd MeshLink/src/mlkernel_geode
make -f Makefile.macosx BUILD=Release

cd MeshLink/src/meshlink
make -f Makefile.macosx BUILD=Release

make -f Makefile machine=macosx BUILD=Release HAVE_GEODE=1 test_harness_c

# Linux

git clone git@github.com:pointwise/MeshLink.git

cd MeshLink/external
tar xzvf [pw-ge]/ge_lite.linux-g++-64.6.8.tgz
ln -s ge_lite.linux-g++-64.6.8/ gelite
cd gelite
ln -s linux-g++-64 linux_x86_64

MeshLink/external/gelite

cd MeshLink/src/mlparser_xerces
make -f Makefile.linux_x86_64 BUILD=Release

cd MeshLink/src/mlkernel_geode
make -f Makefile.linux_x86_64 BUILD=Release

cd MeshLink/src/meshlink
make -f Makefile.linux_x86_64 BUILD=Release

make -f Makefile machine=linux_x86_64 BUILD=Release HAVE_GEODE=1 test_harness_c