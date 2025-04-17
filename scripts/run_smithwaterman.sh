#!/bin/bash

# SmithWaterman build and run script
# This script compiles the C++ implementation of Smith-Waterman and runs it

while getopts "b:r:s:t:o:c:h" opt; do
    case "$opt" in
        b) Build=1 ;;
        r) Run=1 ;;
        s) 
            if [ -z "$seq1" ]; then
                seq1=$OPTARG
            else
                seq2=$OPTARG
            fi
            ;;
        t) threads=$OPTARG ;;
        o) output=$OPTARG ;;
        c) Compare=1 ;;
        h) help=1 ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

# Function to display help message
function show_usage {
    echo "Usage: ./run_smithwaterman.sh [-b] [-r] [-s <seq1>] [-s <seq2>] [-t <threads>] [-o <output>] [-c] [-h]"
    echo "  -b        : Compile the program"
    echo "  -r        : Run the program with the provided parameters"
    echo "  -s        : Path to the sequence files (must provide twice)"
    echo "  -t        : Number of threads for execution"
    echo "  -o        : Output file for results (default: result.txt)"
    echo "  -c        : Compare performance between current and previous runs"
    echo "  -h        : Show this message"
}

if [[ $help -eq 1 ]]; then
    show_usage
    exit 0
fi

# Check for valid arguments
if [[ -z $Build && -z $Run ]]; then
    echo "Error: You must specify at least one of the flags -b or -r"
    show_usage
    exit 1
fi

if [[ $Run -eq 1 && ( -z $seq1 || -z $seq2 ) ]]; then
    echo "Error: When using -r, you must specify -s (twice for both sequence files)"
    show_usage
    exit 1
fi

# Default output file if not provided
if [[ -z $output ]]; then
    output="result.txt"
fi

# Determine the current platform (example only, adjust as needed)
platform="local"

# Path to source file
sourcePath="../src/smithwaterman/smithwaterman.cpp"

# Compiler and flags
compiler="g++"
compileFlags="-std=c++11 -O3 -fopenmp"
outputExe="smithwaterman"

# Compile the program
if [[ $Build -eq 1 ]]; then
    echo "Compiling Smith-Waterman with $compiler..."
    if [[ ! -d "../build" ]]; then
        mkdir ../build
    fi
    compileCommand="$compiler $sourcePath -o ../build/$outputExe $compileFlags"
    echo "Executing: $compileCommand"
    $compileCommand
    if [[ $? -eq 0 ]]; then
        echo "Compilation successful."
    else
        echo "Compilation failed!"
        exit 1
    fi
fi

# Run the program
if [[ $Run -eq 1 ]]; then
    exePath="../build/$outputExe"
    if [[ ! -f $exePath ]]; then
        echo "Error: Executable file does not exist. Please compile first."
        exit 1
    fi
    if [[ ! -f $seq1 || ! -f $seq2 ]]; then
        echo "Error: Sequence files do not exist."
        exit 1
    fi
    runArgs="$seq1 $seq2"
    if [[ $threads -gt 0 ]]; then
        runArgs="$runArgs $threads"
    fi
    runArgs="$runArgs $output"
    echo "Running Smith-Waterman algorithm..."
    $exePath $runArgs
    if [[ $? -eq 0 ]]; then
        echo "Execution successful."
    else
        echo "Error during execution."
        exit 1
    fi
fi
