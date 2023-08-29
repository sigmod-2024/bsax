# bSAX

#  Project Structure

## Code
This folder contains the code for implementing bSAX.

Please refer to file `/sax/include/globals.h` for parameter settings

Place the data in the `dataset` folder and run the following command

```bash
# Move to the code folder
# 1. Compile
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. 
cmake --build .

# 2. Run exp (for example exp5)
./exp5
```

## Data
This folder contains the code for the random dataset generation, as well as the oneDrive links to the real dataset.
These datasets are used in the construction and query answering of bSAX in paper experiments.

OneDrive link: https://1drv.ms/f/s!AoEnTHhjbNSigVaujcVXpN2hqsBF?e=HRzi8r

password: `sigmod2024`


