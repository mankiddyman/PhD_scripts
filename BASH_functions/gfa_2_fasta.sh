gfa_2_fasta() {
    if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
        echo "Usage: gfa_2_fasta input.gfa [output.fasta]"
        return 1
    fi

    input="$1"
    if [ "$#" -eq 2 ]; then
        output="$2"
    else
        # Replace .gfa or .gfa.gz with .fa
        output="${input%.gfa}.fa"
        output="${output%.gfa.gz}.fa"
    fi

    awk '/^S/ {print ">"$2"\n"$3}' "$input" > "$output"
}

