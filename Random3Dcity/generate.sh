
python /Users/ravi/git/Random3Dcity/randomiseCity.py -n 10 -o /Users/ravi/git/mat_util/Random3Dcity/info.xml
python /Users/ravi/git/Random3Dcity/generateCityGML.py -i /Users/ravi/git/mat_util/Random3Dcity/info.xml -o /Users/ravi/git/mat_util/Random3Dcity/CityGML/
python /Users/ravi/git/CityGML2OBJs/CityGML2OBJs.py -i /Users/ravi/git/mat_util/Random3Dcity/CityGML/ -o /Users/ravi/git/mat_util/Random3Dcity/OBJ/
pcl_mesh_sampling /Users/ravi/git/mat_util/Random3Dcity/OBJ/LOD2_0_F0.obj /Users/ravi/git/mat_util/Random3Dcity/PCD/coords.pcd
pcd2npy.py PCD NPY
compute_normals NPY
compute_ma NPY
python /Users/ravi/git/mat_util/region_growing.py NPY