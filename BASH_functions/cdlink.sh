# this function is for when you want to change to the directory of a symbolic link
# it changes to the parent folder of the file.
#
#

#cd for symbolic links
cdlink() {
    if [ -L "$1" ]; then
        cd "$(dirname "$(readlink -f "$1")")" || {
            echo "Error: Failed to change directory." >&2
            return 1
        }
    else
        echo "Error: $1 is not a symbolic link." >&2
        return 1
    fi
}

