//update 1.0.3
pipeline {
agent {label 'merlin'}
    stages {
        stage ("build-local-bwa-flow") {
            steps {
                 dir("ws-bwa-flow") {
                 script {
                    git branch: 'release', url: 'git@github.com:falcon-computing/bwa-flow.git'
                        sh "rm -rf release"
                        sh "mkdir release"
                        dir("release"){
                        sh "rsync -av --exclude=.* /curr/limark/genome-release/build/local/ /curr/limark/falcon2/"
                        sh "rsync -av --exclude=.* /curr/limark/genome-release/build/common/ /curr/limark/falcon2/"
                        sh "source /curr/software/util/modules-tcl/init/bash"
                        sh "module load sdx/17.4; cmake -DCMAKE_BUILD_TYPE=Internal on aws -DRELEASE_VERSION=aws -DDEPLOYMENT_DST=aws -DCMAKE_INSTALL_PREFIX=/curr/limark/falcon2/tools/bin .."
 //                       sh "cmake -DCMAKE_BUILD_TYPE=Release -DRELEASE_VERSION=aws -DDEPLOYMENT_DST= -DCMAKE_INSTALL_PREFIX=/curr/limark/falcon2/tools/bin .."
                        sh "make -j 8"
                        sh "make install"
                        }
                    }
                }
            }
        }
    }
    post {
            always {

                emailext attachLog: true, body: "${currentBuild.currentResult}: Job ${env.JOB_NAME} build ${env.BUILD_NUMBER}\n More info at: ${env.BUILD_URL}console",
                    recipientProviders: [[$class: 'DevelopersRecipientProvider'], [$class: 'RequesterRecipientProvider']],
                    subject: "Jenkins Build ${currentBuild.currentResult}: Job ${env.JOB_NAME}",
                    to: 'udara@limarktech.com, roshantha@limarktech.com'

        }
    }
}
	
