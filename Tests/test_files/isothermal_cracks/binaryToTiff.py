import subprocess
import sys
import time

def install_package(package):
    """Install a package using pip"""
    print(f"Installing {package}...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Try to import required packages, install if missing
try:
    import numpy as np
except ImportError:
    install_package("numpy")
    import numpy as np

try:
    import tifffile as tf
except ImportError:
    install_package("tifffile")
    import tifffile as tf

# Try to import imagecodecs for compression support
compression_available = False
try:
    import imagecodecs
    compression_available = True
    print("Compression support available")
except ImportError:
    try:
        install_package("imagecodecs")
        # Try importing again after installation
        import imagecodecs
        compression_available = True
        print("Compression support installed and available")
    except Exception as e:
        print(f"Could not install or import imagecodecs: {e}")
        print("Will save without compression")
        compression_available = False

print("Starting binary to TIFF conversion...")
start_time = time.time()

# Read dimensions
with open('results/dimensions.txt', 'r') as f:
    dims = f.read().strip().split()
    nX, nY, nZ = int(dims[0]), int(dims[1]), int(dims[2])

print(f"Dimensions: {nX} x {nY} x {nZ}")

# Read binary data
print("Loading binary data...")
binary_file = f'results/microstructure_nX{nX}_nY{nY}_nZ{nZ}.bin'
image_array = np.fromfile(binary_file, dtype=np.uint16)

# Reshape to 3D
image_array = image_array.reshape((nZ, nY, nX))

# Save as TIFF
print("Saving TIFF...")
tiff_file = f'results/microstructure_nX{nX}_nY{nY}_nZ{nZ}.tif'
#tf.imwrite(tiff_file, image_array, compression='lzw')

if compression_available:
    try:
        tf.imwrite(tiff_file, image_array, compression='lzw')
        print("Saved with LZW compression")
    except Exception as e:
        print(f"Compression failed: {e}")
        print("Saving without compression...")
        tf.imwrite(tiff_file, image_array)
else:
    tf.imwrite(tiff_file, image_array)
    print("Saved without compression")

end_time = time.time()
print(f"TIFF conversion completed in {end_time - start_time:.2f} seconds")
print(f"TIFF file saved: {tiff_file}")