#!/bin/sh
CLASSPATH=src:test:config:data

for f in lib/*.jar; do
    CLASSPATH=$CLASSPATH:$f
done

java -Xmx1G -cp $CLASSPATH jline.ConsoleRunner clojure.main -r
