"""The nbconverted rst files have .%5C in the filenames, which is not a valid for
sphinx rst. This script fixes the filenames and removes the .%5C from the
filenames.
"""

from pathlib import Path

# change directory where nbconvert output if changes
directory = Path(".")

for filename in directory.glob("**/*.rst"):
    with open(filename, "r", encoding='utf-8') as file:
        content = file.read()

    if ".%5C" in content:
        new_content = content.replace(".%5C", "")

        with open(filename, "w", encoding="utf-8") as file:
            file.write(new_content)