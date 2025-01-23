#!/bin/bash

# Define the file path
file_path="nmmo_x_0_3.data"

# Define the text to be replaced and the replacement text
text_to_replace="59       atoms\n3  atom types"
replacement_text="59 atoms\n0 bonds\n0 angles\n4 atom types\n0 bond types\n0 angle types"

# Use sed to replace the text
sed -i "s/${text_to_replace}/${replacement_text}/g" "$file_path"

# Define the text to be added before 'Atoms'
text_to_add="59 atoms\n0 bonds\n0 angles\n4 atom types\n0 bond types\n0 angle types"

# Use sed to add text before 'Atoms'
sed -i "/Atoms/ i ${text_to_add}" "$file_path"

echo "Text replacement and addition complete."

