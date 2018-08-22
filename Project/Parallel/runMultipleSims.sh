echo "Running multiple simulations..."

make main



bash cleanData.sh
./main nx 16 ny 16 endTime 0 sigma 10000000000
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/1616/Initial/

bash cleanData.sh
./main nx 16 ny 16 endTime 10 sigma 10000000000
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/1616/

bash cleanData.sh
./main nx 32 ny 32 endTime 0 sigma 10000000000
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/3232/Initial/

bash cleanData.sh
./main nx 32 ny 32 endTime 10 sigma 10000000000
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/3232/

bash cleanData.sh
./main nx 64 ny 64 endTime 0 sigma 10000000000
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/6464/Initial/

bash cleanData.sh
./main nx 64 ny 64 endTime 10 sigma 10000000000
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/6464/

bash cleanData.sh
./main nx 128 ny 128 endTime 0 sigma 10000000000
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/128128/Initial/

bash cleanData.sh
./main nx 128 ny 128 endTime 10 sigma 10000000000
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/128128/

bash cleanData.sh
./main nx 256 ny 256 endTime 0 sigma 10000000000
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/256256/Initial/

bash cleanData.sh
./main nx 256 ny 256 endTime 10 sigma 10000000000
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/256256/

bash cleanData.sh
./main nx 512 ny 512 endTime 0 sigma 10000000000
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/512512/Initial/

bash cleanData.sh
./main nx 512 ny 512 endTime 10 sigma 10000000000
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionIdealLimit/512512/

bash cleanData.sh
./main nx 16 ny 16 endTime 0 sigma 100
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/1616/Initial/

bash cleanData.sh
./main nx 16 ny 16 endTime 10 sigma 100
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/1616/

bash cleanData.sh
./main nx 32 ny 32 endTime 0 sigma 100
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/3232/Initial/

bash cleanData.sh
./main nx 32 ny 32 endTime 10 sigma 100
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/3232/

bash cleanData.sh
./main nx 64 ny 64 endTime 0 sigma 100
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/6464/Initial/

bash cleanData.sh
./main nx 64 ny 64 endTime 10 sigma 100
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/6464/

bash cleanData.sh
./main nx 128 ny 128 endTime 0 sigma 100
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/128128/Initial/

bash cleanData.sh
./main nx 128 ny 128 endTime 10 sigma 100
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/128128/

bash cleanData.sh
./main nx 256 ny 256 endTime 0 sigma 100
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/256256/Initial/

bash cleanData.sh
./main nx 256 ny 256 endTime 10 sigma 100
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/256256/

bash cleanData.sh
./main nx 512 ny 512 endTime 0 sigma 100
cp -r Data/Final/* ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/512512/Initial/

bash cleanData.sh
./main nx 512 ny 512 endTime 10 sigma 100
cp -r Data/Final ./../../../DataToKeep/CUDAPaper/Convergence/FieldLoopAdvectionSigma100/512512/
