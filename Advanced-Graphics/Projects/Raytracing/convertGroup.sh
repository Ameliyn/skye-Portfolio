#!/bin/bash

for name in texturemovie/*;
do foo="${name%.xwd}".jpg;
   foo=${foo#"texturemovie/"};
   foo=texturemoviejpg/"$foo";
   convert $name $foo ;
done
