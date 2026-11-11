def toolVersion(command: String) {
    return "{ ${command}; } 2>/dev/null | grep -vi 'not found' | grep -m1 . || echo ''"
}

def rVersion() {
    return toolVersion("R --version | head -n1 | sed -E 's/R version ([0-9.]+).*/\\1/'")
}

def rPackageVersion(package_name: String) {
    return toolVersion("Rscript -e 'cat(as.character(packageVersion(commandArgs(TRUE)[1])))' ${package_name}")
}

def perlVersion() {
    return toolVersion("perl --version | grep -Eo 'v[0-9.]+' | head -n1 | sed 's/v//'")
}

def samtoolsVersion() {
    return toolVersion("samtools --version | head -n1 | sed 's/samtools //'")
}

def pigzVersion() {
    return toolVersion("pigz --version 2>&1 | sed 's/pigz //'")
}
