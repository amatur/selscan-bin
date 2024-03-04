
To install: clone the repository

`git clone https://github.com/amatur/selscan-bin`     

Then, run the makefile

From linux
`make     `

From Mac/OSX
`make -f Makefile.osx   `


To compute iHS
` ./selbin ehh --all --vcf example/p4.2000.vcf`     

To compute EHH (Here loc is the order of the position, not actual physical position)
` ./selbin ehh --loc 0 --vcf example/p4.2000.vcf`     


