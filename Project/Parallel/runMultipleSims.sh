echo "Running multiple simulations..."

make main



################################################################################
#########  gamma 4/3 with By=0.1
################################################################################

gammanum=4
gammaden=3


bash cleanData.sh
./main nx 128 ny 128 frameSkip 12 gammanum $gammanum gammaden=$gammaden
cp -r Data/* ./../../../DataToKeep/KHI/Gam5by3B0p1/Data128128/

bash cleanData.sh
./main nx 256 ny 256 frameSkip 48 gammanum $gammanum gammaden=$gammaden
cp -r Data/* ./../../../DataToKeep/KHI/Gam5by3B0p1/Data256256/

bash cleanData.sh
./main nx 512 ny 512 frameSkip 192 gammanum $gammanum gammaden=$gammaden
cp -r Data/* ./../../../DataToKeep/KHI/Gam5by3B0p1/Data512512/

echo "Finished!"
