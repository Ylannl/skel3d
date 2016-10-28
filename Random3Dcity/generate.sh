mkdir CityGML
mkdir OBJ
mkdir PCD
mkdir NPY
python /Users/ravi/git/Random3Dcity/randomiseCity.py -p 1 -n 16 -o /Users/ravi/git/mat_util/Random3Dcity/info.xml
python /Users/ravi/git/Random3Dcity/generateCityGML.py -i /Users/ravi/git/mat_util/Random3Dcity/info.xml -o /Users/ravi/git/mat_util/Random3Dcity/CityGML/
python /Users/ravi/git/CityGML2OBJs/CityGML2OBJs.py -i /Users/ravi/git/mat_util/Random3Dcity/CityGML/ -o /Users/ravi/git/mat_util/Random3Dcity/OBJ/
pcl_mesh_sampling -leaf_size 0.1 -n_samples 200000 /Users/ravi/git/mat_util/Random3Dcity/OBJ/LOD2_1_F0_DSM.obj /Users/ravi/git/mat_util/Random3Dcity/PCD/coords.pcd
python /Users/ravi/git/mat_util/scripts/io/pcd2npy.py PCD NPY
compute_normals NPY
compute_ma NPY
python3 /Users/ravi/git/mat_util/scripts/segment.py NPY
