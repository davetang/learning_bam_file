## README

If you are using macOS and want to build the README.md using the Makefile, you can use Docker and this [Ubuntu image](https://hub.docker.com/r/davetang/build).

```bash
docker pull davetang/build:1.0
```

However, when I mounted the learning_bam_file directory to the Docker container and ran the Makefile from there, I'd get an [eterm](https://github.com/conda/conda/issues/6603) error. (Specifically, I got an error with ncurses-6.2-he6710b0_1/share/terminfo/E/Eterm-color.) It turns out that for macOS High Sierra (or later), the default file system is APFS and [it is case-insensitive by default](https://docker-docs.netlify.app/docker-for-mac/osxfs/) and Docker inherits this! Therefore, when you use Docker to run the Makefile, don't run it in a mounted volumne. For example:

```bash
docker run --rm -it -v $(pwd):/work davetang/build:1.0 /bin/bash

# inside the Docker container
cd /tmp
git clone https://github.com/davetang/learning_bam_file.git
cd learning_bam_file
make
```

