## Assignment for Week 1 ##

Practicing creating and adding a file to my repository from the terminal:

```bash
echo "Hello world" > assignment1.txt
git add assignment1.txt
git commit -m "Add assignment1.txt"
git push origin main
```

Checking samtools version for bioinfo env:
```bash
cd ~/bioinfo
pwd
samtools --version
samtools 1.22.1
Using htslib 1.22.1
Copyright (C) 2025 Genome Research Ltd.
```

Creating and changing to a new directory for week 1
```bash
mkdir assignment_week1
cd ~/bioinfo/assignment_week1
ls
pwd
```
Create nested directories inside assignment_week1 directory:

```bash
mkdir homework1
cd homework1
pwd
mkdir data 
cd data
pwd
```
Create nested directories at once with a single command: 
```bash
mkdir -p ~/nested1/nested2/nested3
```
Change to the root directory using an absolute path:
```bash
cd /Users/jessicabuchser
```
vs. using a relative path:

```bash
cd /
cd /Users
cd /jessicabuchser
```
Navigate to the subdirectory data and then navigate up two directories, then go up two directories at once:

```bash
cd ~/bioinfo/assignment_week1/homework1/data
pwd
cd ..
pwd
cd ../..
pwd
```
Create an empty file in assignment_week1 and move it to a new directory called temp:

```bash
cd ~/bioinfo/assignment_week1
touch newfile.txt
mkdir temp
mv newfile.txt temp
ls temp
```