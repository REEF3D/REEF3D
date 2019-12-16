for i in `find . -name "*.h" -type f`; do
    sed -i '' 's/a->eps/ee/g' $i 
done
