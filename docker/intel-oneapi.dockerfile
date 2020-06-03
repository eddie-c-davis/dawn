FROM intel/oneapi-hpckit:latest

ENV PATH=/opt/intel/inteloneapi/compiler/latest/linux/bin/intel64:/opt/intel/inteloneapi/intelpython/python3.7/bin:/usr/local/bin:$PATH \
    LD_LIBRARY_PATH=/opt/intel/inteloneapi/compiler/latest/linux/lib:/opt/intel/inteloneapi/compiler/latest/linux/compiler/lib/intel64:/opt/intel/inteloneapi/intelpython/python3.7lib:$LD_LIBRARY_PATH \
    mpich_version=3.1.4 \
    OMP_NUM_THREADS=1 \
    BOOST_HOME=/usr/local/boost_1_73_0 CXX=icpc CC=icc CPPFLAGS="-I/usr/local/boost_1_73_0" PYTHON_CFLAGS="-O3 -std=c90"

##
## Install some necessary prerequisities
##

RUN apt-get update && apt-get install -y git wget make cmake curl libssl-dev llvm clang libclang-dev clang-tidy clang-format && \
    cd /usr/local && \
    wget https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.gz && \
    tar xzf boost_1_73_0.tar.gz && \
    rm boost_1_73_0.tar.gz 

##
## Build MPICH. Update environment to use MPICH instead of Intel MPI
##

#RUN wget -q http://www.mpich.org/static/downloads/${mpich_version}/mpich-${mpich_version}.tar.gz && \
#    tar xzf mpich-${mpich_version}.tar.gz && \
#    cd mpich-${mpich_version} && \
#    CC=icc CXX=icpc FC=ifort ./configure --prefix=/usr/local --enable-cxx --enable-fortran --enable-fast=all,O3 && \
#    make -j4 && make install && ldconfig && \
#    cd .. && rm -rf mpich-${mpich_version}*
#
#ENV PATH=/usr/local/bin:$PATH \
#    LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH \

##
## Install mpi4py with MPICH built previously
##

#RUN python -m pip uninstall mpi4py && \
#    MPICC=/usr/local/bin/mpicc python -m pip install mpi4py
     
# ---------------------- GridTools ----------------------
# RUN curl -L https://github.com/GridTools/gridtools/archive/v1.0.4.tar.gz | \
#     tar -xz -C /usr/src && \
#     mkdir /usr/src/gridtools-1.0.4/build
# RUN cmake -S /usr/src/gridtools-1.0.4 -B /usr/src/gridtools-1.0.4/build \
#     -DBUILD_TESTING=OFF \
#     -DINSTALL_TOOLS=OFF \
#     -DGT_INSTALL_EXAMPLES=OFF \
#     -DCMAKE_BUILD_TYPE=Release \
#     -DCMAKE_INSTALL_PREFIX=/usr/local/gridtools \
#     -GNinja && \
#     cmake --build /usr/src/gridtools-1.0.4/build -j $(nproc) --target install && \
#     rm -rf /usr/src/gridtools-1.0.4/build
# Other python dependencies for using and testing dawn
RUN python -m pip install pytest
RUN python -m pip install --no-cache-dir scikit-build
# CMake
RUN apt purge -y cmake
RUN cd /tmp && \
    wget https://cmake.org/files/v3.16/cmake-3.16.5.tar.gz && \
    tar -xzvf cmake-3.16.5.tar.gz && \
    cd cmake-3.16.5/ && \
    ./bootstrap && \
    make -j$(nproc) && \
    make install && \
    rm -rf cmake-3.16.5/ && \
    cd /
# ---------------------- Protobuf ----------------------
# RUN curl -L https://github.com/protocolbuffers/protobuf/releases/download/v3.10.1/protobuf-all-3.10.1.tar.gz | \
#     tar -xz -C /usr/src
# # These files seem to have a high UID/GID by default, so update this
# RUN chown root:root /usr/src/protobuf-3.10.1 -R
# RUN cmake -S /usr/src/protobuf-3.10.1/cmake -B /usr/src/protobuf-3.10.1/build \
#     -Dprotobuf_BUILD_EXAMPLES=OFF \
#     -Dprotobuf_BUILD_TESTS=OFF \
#     -Dprotobuf_INSTALL_EXAMPLES=OFF \
#     -Dprotobuf_BUILD_PROTOC_BINARIES=ON \
#     -DCMAKE_BUILD_TYPE=Release \
#     -DCMAKE_INSTALL_PREFIX=/usr/local/protobuf \
#     -DBUILD_SHARED_LIBS=ON \
#     -GNinja && \
#     cmake --build /usr/src/protobuf-3.10.1/build --target install -j $(nproc) && \
#     rm -rf /usr/src/protobuf-3.10.1/build
# RUN cd /usr/src/protobuf-3.10.1/python && \
#     PROTOC=/usr/local/protobuf/bin/protoc python setup.py build && \
#     mkdir -p /usr/local/protobuf/lib/python && \
#     mv /usr/src/protobuf-3.10.1/python/build/lib/google /usr/local/protobuf/lib/python/google
# ---------------------- Dawn ----------------------
RUN git clone https://github.com/MeteoSwiss-APN/dawn.git /usr/src/dawn && \ 
    CC=icc CXX=icpc && \
    mkdir /usr/src/dawn/build && \
    cd /usr/src/dawn/build && \
    cmake .. -B . -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local && \
    cmake --build . --target install -j $(nproc)

##
## Install GT4Py
##
# RUN git clone --depth 1 -b fv3_validation https://github.com/eddie-c-davis/gt4py.git /usr/src/gt4py && \
#     CC=icc CXX=icpc CPPFLAGS="-I/usr/local/boost_1_73_0" python -m pip install --no-cache-dir --no-binary ":all:" -e /usr/src/gt4py && \
#     cd /usr/src/gt4py && python -m gt4py.gt_src_manager install
