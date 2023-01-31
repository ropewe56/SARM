curl https://sh.rustup.rs -sSf > rustup.sh
bash  rustup.sh

rustc --version
rustup completions bash > /etc/bash_completion.d/rustup.bash-completion
rustup update
rustup self update

rustup install nightly
rustup default nightly