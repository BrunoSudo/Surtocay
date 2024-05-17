#!/bin/bash

#!/bin/bash

# Specify the source files
SOURCE_FILES=($(find . -name "*.cpp" -not -name "tester.cpp") "../eos_base.cpp" "../eos_hotQCD.cpp"  "../pretty_ostream.cpp" "../util.cpp" "../emoji.cpp")

# Compile the program
if g++ -o "Surtocay" "${SOURCE_FILES[@]}" -lgsl -lgslcblas -lm; then
    echo "Compilation successful."

    # Run the program if compilation was successful
    if [ -f "Surtocay" ]; then
        echo "Running Surtocay..."
        ./Surtocay
    else
        echo "Executable not found."
    fi
else
    echo "Compilation failed."
fi


