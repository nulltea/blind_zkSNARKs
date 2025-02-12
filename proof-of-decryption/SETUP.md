# Setup
This setup is only intended for hosts with Intel or x86 architecture. It is assumed that you have a decent C compiler installed on your machine.

## Prerequisite  
Install these tools.  

```bash
sudo apt-get update
sudo apt-get install libmpfr-dev
sudo apt-get install libgmp-dev
pip3 install cffi
```

#### Availability of AVX512 
This setup is performant if your CPU supports AVX512. To check this, you can run the following command.
If you find AVX512 or any mention of it in the output of the second command in the `Flags` section, your CPU supports it.  

```bash
lscpu
grep -o 'avx512[^ ]*' /proc/cpuinfo
```

### Installation
On the root folder of the project, execute:  

```bash
# For a basic install
make

# For an installation with AVX512 support
make all
```


For the **Python** project:    

```bash
cd python
make
```

