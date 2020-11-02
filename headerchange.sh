for i in `find . -name "*.cpp" -type f`; do
    sed -i '' 's/Copyright 2008-2020 Hans Bihs/Copyright 2018-2020 Tobias Martin/g' $i 
done
