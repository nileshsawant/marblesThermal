import subprocess
import sys
import time
import os

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

try:
    from skimage import measure
except ImportError:
    install_package("scikit-image")
    from skimage import measure

try:
    import trimesh
except ImportError:
    install_package("trimesh")
    import trimesh

# Optional: For STEP/IGES export (requires additional dependencies)
try:
    import FreeCAD
    import Part
    freecad_available = True
    print("FreeCAD available for STEP/IGES export")
except ImportError:
    freecad_available = False
    print("FreeCAD not available - only STL export supported")

def find_tiff_file():
    """Find the TIFF file in results directory or current directory"""
    # Check current directory first
    current_dir = "."
    tiff_files = [f for f in os.listdir(current_dir) if f.endswith('.tif')]
    
    if tiff_files:
        tiff_file = tiff_files[0]
        if len(tiff_files) > 1:
            print(f"Multiple TIFF files found in current directory: {tiff_files}")
            print(f"Using: {tiff_file}")
        return tiff_file
    
    # Check results directory
    results_dir = "results"
    if os.path.exists(results_dir):
        tiff_files = [f for f in os.listdir(results_dir) if f.endswith('.tif')]
        if tiff_files:
            tiff_file = os.path.join(results_dir, tiff_files[0])
            if len(tiff_files) > 1:
                print(f"Multiple TIFF files found in results directory: {tiff_files}")
                print(f"Using: {tiff_files[0]}")
            return tiff_file
    
    raise FileNotFoundError("No TIFF files found in current directory or results/")

def generate_mesh_for_tag(volume_data, tag_value, smoothing=True):
    """Generate mesh for a specific tag value"""
    if tag_value == -1:
        # Special case for combined mesh
        binary_mask = (volume_data > 0).astype(np.uint8)
        print(f"Generating combined mesh for all solid materials...")
    else:
        print(f"Generating mesh for tag {tag_value}...")
        binary_mask = (volume_data == tag_value).astype(np.uint8)
    
    # Check if tag exists
    if np.sum(binary_mask) == 0:
        print(f"Warning: No voxels found for tag {tag_value if tag_value != -1 else 'combined'}")
        return None
    
    print(f"Found {np.sum(binary_mask)} voxels for tag {tag_value if tag_value != -1 else 'combined'}")
    
    # Generate mesh using marching cubes
    try:
        verts, faces, normals, values = measure.marching_cubes(
            binary_mask, 
            level=0.5,
            spacing=(1.0, 1.0, 1.0)
        )
        
        print(f"Generated mesh: {len(verts)} vertices, {len(faces)} faces")
        
        # Create trimesh object
        mesh = trimesh.Trimesh(vertices=verts, faces=faces, vertex_normals=normals)
        
        # Apply smoothing using correct trimesh API
        if smoothing and len(verts) > 0:
            try:
                # Use smoothing filter if available
                mesh = mesh.smoothed()
            except AttributeError:
                try:
                    # Alternative smoothing method
                    mesh = trimesh.smoothing.filter_laplacian(mesh, iterations=1)
                except (AttributeError, ImportError):
                    # Skip smoothing if not available
                    print("Smoothing not available - using original mesh")
        
        # Clean up mesh
        mesh.remove_degenerate_faces()
        mesh.remove_duplicate_faces()
        mesh.remove_unreferenced_vertices()
        
        if not mesh.is_empty:
            print(f"Final mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
            return mesh
        else:
            print(f"Warning: Empty mesh generated for tag {tag_value}")
            return None
            
    except Exception as e:
        print(f"Error generating mesh for tag {tag_value}: {e}")
        return None

def export_stl(mesh, filename, tag_value):
    """Export mesh to STL format"""
    try:
        if tag_value == -1:
            stl_filename = f"{filename}_combined.stl"
        else:
            stl_filename = f"{filename}_tag{tag_value}.stl"
        mesh.export(stl_filename)
        print(f"Exported STL: {stl_filename}")
        return stl_filename
    except Exception as e:
        print(f"Error exporting STL for tag {tag_value}: {e}")
        return None

def export_step_iges_freecad(mesh, filename, tag_value, format='step'):
    """Export mesh to STEP/IGES using FreeCAD"""
    if not freecad_available:
        print(f"FreeCAD not available - cannot export {format.upper()}")
        return None
    
    try:
        # Create FreeCAD document
        doc = FreeCAD.newDocument()
        
        # Convert trimesh to FreeCAD mesh
        import Mesh
        freecad_mesh = Mesh.Mesh()
        
        # Add triangles to FreeCAD mesh
        for face in mesh.faces:
            triangle = [mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]]
            freecad_mesh.addFacet(triangle[0], triangle[1], triangle[2])
        
        # Create mesh object in document
        mesh_obj = doc.addObject("Mesh::Feature", f"Tag{tag_value}")
        mesh_obj.Mesh = freecad_mesh
        
        # Export based on format
        if format.lower() == 'step':
            if tag_value == -1:
                output_filename = f"{filename}_combined.step"
            else:
                output_filename = f"{filename}_tag{tag_value}.step"
                
        elif format.lower() == 'iges':
            if tag_value == -1:
                output_filename = f"{filename}_combined.iges"
            else:
                output_filename = f"{filename}_tag{tag_value}.iges"
            Part.export([mesh_obj], output_filename)
        
        # Close document
        FreeCAD.closeDocument(doc.Name)
        
        print(f"Exported {format.upper()}: {output_filename}")
        return output_filename
        
    except Exception as e:
        print(f"Error exporting {format.upper()} for tag {tag_value}: {e}")
        return None

def main():
    print("Starting TIFF to CAD conversion...")
    start_time = time.time()
    
    try:
        # Find TIFF file
        tiff_file = find_tiff_file()
        print(f"Loading TIFF file: {tiff_file}")
        
        # Load TIFF data
        volume_data = tf.imread(tiff_file)
        print(f"TIFF dimensions: {volume_data.shape}")
        print(f"Data type: {volume_data.dtype}")
        
        # Get unique tags (excluding 0 for void)
        unique_tags = np.unique(volume_data)
        unique_tags = unique_tags[unique_tags > 0]  # Remove void (0)
        print(f"Found tags: {unique_tags}")
        
        if len(unique_tags) == 0:
            print("No non-zero tags found in TIFF file")
            return
        
        # Create output directory
        output_dir = "cad_exports"
        os.makedirs(output_dir, exist_ok=True)
        
        # Base filename
        base_filename = os.path.splitext(os.path.basename(tiff_file))[0]
        
        exported_files = []
        
        # Process each tag (limit to first few tags for testing)
        max_tags_to_process = min(100, len(unique_tags))  # Process max 100 tags for testing
        print(f"Processing first {max_tags_to_process} tags (out of {len(unique_tags)}) for faster testing...")
        
        for i, tag in enumerate(unique_tags[:max_tags_to_process]):
            print(f"\nProcessing tag {tag} ({i+1}/{max_tags_to_process})...")
            
            # Generate mesh
            mesh = generate_mesh_for_tag(volume_data, tag, smoothing=False)  # Disable smoothing for speed
            
            if mesh is not None:
                output_base = os.path.join(output_dir, base_filename)
                
                # Export STL (always available)
                stl_file = export_stl(mesh, output_base, tag)
                if stl_file:
                    exported_files.append(stl_file)
                
                # Export STEP/IGES if FreeCAD is available
                if freecad_available:
                    step_file = export_step_iges_freecad(mesh, output_base, tag, 'step')
                    if step_file:
                        exported_files.append(step_file)
                    
                    iges_file = export_step_iges_freecad(mesh, output_base, tag, 'iges')
                    if iges_file:
                        exported_files.append(iges_file)
                else:
                    print("For STEP/IGES export, install FreeCAD:")
                    print("conda install -c conda-forge freecad")
        
        # Export combined mesh (all tags together)
        print(f"\nGenerating combined mesh...")
        combined_mesh = generate_mesh_for_tag(volume_data, -1, smoothing=False)  # Special case for combined
        
        if combined_mesh is not None:
            output_base = os.path.join(output_dir, base_filename)
            
            # Export combined STL
            stl_file = export_stl(combined_mesh, output_base, -1)
            if stl_file:
                exported_files.append(stl_file)
            
            # Export combined STEP/IGES if available
            if freecad_available:
                step_file = export_step_iges_freecad(combined_mesh, output_base, -1, 'step')
                if step_file:
                    exported_files.append(step_file)
        
        # Summary
        end_time = time.time()
        print(f"\nCAD conversion completed in {end_time - start_time:.2f} seconds")
        print(f"Exported {len(exported_files)} files:")
        for file in exported_files:
            if os.path.exists(file):
                file_size = os.path.getsize(file) / (1024*1024)  # MB
                print(f"  - {file} ({file_size:.1f} MB)")
            else:
                print(f"  - {file} (file not found)")
        
        print(f"\nFiles saved in: {output_dir}/")
        print("\nCAD files can be opened in:")
        print("- STL: Blender, MeshLab, FreeCAD, SolidWorks, Fusion 360")
        print("- STEP: FreeCAD, SolidWorks, Fusion 360, Onshape")
        print("- IGES: FreeCAD, SolidWorks, Fusion 360")
        
        if len(unique_tags) > max_tags_to_process:
            print(f"\nNote: Only processed {max_tags_to_process} out of {len(unique_tags)} tags.")
            print("To process all tags, change max_tags_to_process in the script.")
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()