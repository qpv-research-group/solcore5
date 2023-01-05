echo "Installing gfortran for arm64"
curl -L https://github.com/fxcoudert/gfortran-for-macOS/releases/download/12.1-monterey/gfortran-ARM-12.1-Monterey.dmg -o gfortran.dmg
GFORTRAN_SHA256=$(shasum -a 256 gfortran.dmg)
KNOWN_SHA256="e2e32f491303a00092921baebac7ffb7ae98de4ca82ebbe9e6a866dd8501acdf  gfortran.dmg"

if [ "$GFORTRAN_SHA256" != "$KNOWN_SHA256" ]; then
    echo sha256 mismatch
    exit 1
fi

hdiutil attach -mountpoint /Volumes/gfortran gfortran.dmg
sudo installer -pkg /Volumes/gfortran/gfortran.pkg -target /
sudo mv /usr/local/bin/gfortran /usr/local/bin/gfortran_arm64
type -p gfortran_arm64


echo "Installing gfortran for x86_64"
brew reinstall gfortran
type -p gfortran

echo '\n Installing ninja'
brew install ninja

echo '\n Installing build system'
pip install meson-python cython numpy