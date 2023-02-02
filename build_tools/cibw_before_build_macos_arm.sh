#!/bin/bash

echo "++++++ Installing arm64_cross_gfortran"
curl -L -O https://github.com/isuruf/gcc/releases/download/gcc-11.3.0-2/gfortran-darwin-arm64-cross.tar.gz
export GFORTRAN_SHA=527232845abc5af21f21ceacc46fb19c190fe804
if [[ "$(shasum gfortran-darwin-arm64-cross.tar.gz)" != "${GFORTRAN_SHA}  gfortran-darwin-arm64-cross.tar.gz" ]]; then
    echo "shasum mismatch for gfortran-darwin-arm64-cross"
    exit 1
fi
sudo mkdir -p /opt/
sudo cp "gfortran-darwin-arm64-cross.tar.gz" /opt/gfortran-darwin-arm64-cross.tar.gz
pushd /opt
    sudo tar -xvf gfortran-darwin-arm64-cross.tar.gz
    sudo rm gfortran-darwin-arm64-cross.tar.gz
popd
export FC_ARM64="$(find /opt/gfortran-darwin-arm64-cross/bin -name "*-gfortran")"

export FC_LOC=/opt/gfortran-darwin-arm64-cross/bin
libgfortran="$(find /opt/gfortran-darwin-arm64-cross/lib -name libgfortran.dylib)"
libdir="$(dirname "$libgfortran")"

export FC_LIBDIR=$libdir
export FC_ARM64_LDFLAGS="-L$libdir -Wl,-rpath,$libdir"
if [[ "${PLAT:-}" == "arm64" ]]; then
    export FC=$FC_ARM64
fi

which gfortran
type -p gfortran

/opt/gfortran-darwin-arm64-cross/bin/arm64-apple-darwin20.0.0-gfortran --version

sudo xcode-select -switch /Applications/Xcode_12.5.1.app
export SDKROOT=/Applications/Xcode_12.5.1.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX11.3.sdk

printf '\n++++++ Re-installing gfortran x86_64 \n'
brew reinstall gfortran
type -p gfortran
type -p gcc

printf '\n++++++ Installing ninja \n '
brew install ninja

printf '\n++++++ Installing build system \n'
pip install meson-python cython numpy