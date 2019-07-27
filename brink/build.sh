lib_root=$HOME/documents/processing/libraries/brink
files=`find -f $lib_root/src \( -name "*.java" \)`

cd $lib_root/src
# for file in $files ; do
    javac -d . -classpath ~/applications/processing.app/contents/java/core/library/core.jar $files
# done
jar cvf ../library/brink.jar brink
find -f $lib_root/src \( -name "*.class" \) | xargs rm -f
