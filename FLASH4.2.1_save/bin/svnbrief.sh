#!/bin/sh -f

# script to create a brief summary of svn situation
# to append to FLASH release


rawversions=$(echo \
 `svn status -v ../source ../bin|grep -v '^[?]  '|grep -v '^Status '|cut -b 9-|sed -e 's/^ *//'| \
  sed -e '/-/d' -e 's/ .*//'|sort -n -r|sed -e 's/^/r/'|uniq -c|sort -n -r|sed -e '1s/^ *[0-9]* */svn:/'| \
  tr -d ' '`|tr ' ' ',')

versions=$rawversions

info=$(echo `svn st ../source ../bin|grep -v '^[?]  '|cut -b 1|sort|uniq -c|tr -d ' '`|tr ' ' ',')
if [ -n "$info" ]; then
    versions="$rawversions,changed:$info"
fi

if [ -n "$versions" ]; then
    echo "($versions)"
fi
