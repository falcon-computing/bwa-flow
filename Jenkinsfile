//update 1.0.9
pipeline {
agent {label 'merlin'}
    stages {
        stage ("build-aws-bwa-flow") {
            steps {
                 dir("ws-bwa-flow") {
		             checkout scm
                 script {
                        sh "rm -rf release"
                        sh "mkdir release"
                        dir("release"){
                        sh "source /curr/software/util/modules-tcl/init/bash"
                        version= sh(returnStdout: true, script: 'git describe --tag').trim()
                        sh "echo $version"
                        sh "module load sdx/17.4; cmake -DCMAKE_BUILD_TYPE=Release -DRELEASE_VERSION=$version -DDEPLOYMENT_DST=aws -DCMAKE_INSTALL_PREFIX=/curr/limark/falcon2/tools/bin .."
                        sh "make -j 8"
			sh "make test"
                        sh "make install"
                        sh "cd ~/falcon2/tools/bin; echo s3://fcs-cicd-test/release/aws/bwa-flow/bwa-flow-$version-aws > latest"
                        sh "cd ~/falcon2/tools/bin; aws s3 cp bwa-flow s3://fcs-cicd-test/release/aws/bwa-flow/bwa-flow-$version-aws"
                        sh "cd ~/falcon2/tools/bin; aws s3 cp latest s3://fcs-cicd-test/release/aws/bwa-flow/latest"
                        sh "cd ~/falcon2/tools/bin; rm -f latest"
                        }
                    }
                }
            }
        }
    }
}
