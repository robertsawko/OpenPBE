sed -i -e 's/binary/ascii/' -e 's/ compressed/ uncompressed/' system/controlDict
foamFormatConvert
blockMesh
topoSet
subsetMesh -overwrite half
sed -i -e 's/oldInternalFaces/symmetry1/' 0/* constant/polyMesh/boundary
topoSet
subsetMesh -overwrite quarter
sed -i -e 's/oldInternalFaces/symmetry2/' 0/* constant/polyMesh/boundary
sed -i -e 's/empty/symmetryPlane/' 0/* constant/polyMesh/boundary
sed -i -e 's/ascii/binary/' -e 's/ uncompressed/ compressed/' system/controlDict
foamFormatConvert