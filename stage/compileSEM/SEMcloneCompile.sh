git clone https://github.com/sem3d/SEM.git
cd SEM
mkdir buildSEM
mkdir buildRF
cd ..
chmod 777 compileSEMAmu.sh
chmod 777 compileRFAmu.sh
mv compileSEMAmu.sh SEM/buildSEM
mv compileRFAmu.sh SEM/buildRF
mv makefileRFamu SEM/buildRF
cd SEM/buildRF
./compileRFAmu.sh
cd ../buildSEM
./compileSEMAmu.sh
cd ../../
