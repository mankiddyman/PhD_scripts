#set a fancy prompt (non-color, unless we know we "want" color)

case "$TERM" in
	    xterm-color|*-256color) color_prompt=yes;;
	        xterm) color_prompt=yes;;
	esac
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'
    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

#getting modules working
source /opt/share/software/scs/appStore/modules/init/profile.sh



#adding favourite directories
netscratch_home=~/../../netscratch/dep_mercier/grp_marques/Aaryan

export PATH=$PATH:~/julia-1.9.3/bin
export PATH=$PATH:~/nextflow

# for copy pasting
#to load modules
source /netscratch/common/Saurabh/envMod/init/bash
#colours terminal prompt
export PS1="\[\e[32m\]\u@\h\[\e[0m\]:\[\e[34m\]\w\[\e[0m\] \[\e[91m\]\t\[\e[0m\]\$ "
#colour ls of sh files differently
export LS_COLORS="$LS_COLORS:*.sh=01;35"
export PATH=$PATH:/netscratch/dep_mercier/grp_marques/Aaryan/gffread/
export PATH=$PATH:/netscratch/dep_mercier/grp_marques/Aaryan/methods/HapHiC/
export PATH=$PATH:/netscratch/dep_mercier/grp_marques/Aaryan/methods/bioawk/
export PATH=$PATH:/netscratch/dep_mercier/grp_marques/Aaryan/methods/KMC/bin
export PATH=$HOME/.local/bin:$PATH
export PATH=$HOME:/netscratch/dep_mercier/grp_marques/Aaryan/methods/fastoche/bin:$PATH
export PATH=$HOME:/netscratch/dep_mercier/grp_marques/Aaryan/methods/bin:$PATH
export MERQURY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/merqury/share/merqury




# >>> mamba initialize >>>
# !! Contents within this block are managed by 'micromamba shell init' !!
export MAMBA_EXE='/home/abhatia/.local/bin/micromamba';
export MAMBA_ROOT_PREFIX='/home/abhatia/micromamba';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias micromamba="$MAMBA_EXE"  # Fallback on help from micromamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<

#new command for fasta summarising , prints chromsome names and lengths
fasta_summary() {
	    if [[ "$1" == *.gz ]]; then
		            zcat "$1" | awk '/^>/ {if (seq_length) print header, seq_length; header=$0; seq_length=0; next} {seq_length+=length($0)} END {if (seq_length) print header, seq_length}'
			        else
					        cat "$1" | awk '/^>/ {if (seq_length) print header, seq_length; header=$0; seq_length=0; next} {seq_length+=length($0)} END {if (seq_length) print header, seq_length}'
						    fi
					    }


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

export PATH="$PATH:~/Downloads/vivaldi/opt/vivaldi"
export PATH="$PATH:~/Downloads/sublime_text/opt/sublime_text" 
export PATH="$PATH:~/Downloads"
export PATH=$HOME/local/bin:$PATH
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"  # This loads nvm
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"  # This loads nvm bash_completion
. "$HOME/.cargo/env"

alias nvim="$HOME/Downloads/nvim_folder/usr/bin/nvim"


export TERM=xterm-256color
