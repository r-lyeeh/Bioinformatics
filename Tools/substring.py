import re
f= open("input.txt")
count = 0
for index in f:

str1=""
str2=""
a=[m.start() for m in re.finditer(str2,str1)]
b=[i+1 for i in a]
print b
