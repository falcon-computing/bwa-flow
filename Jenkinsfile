//update 1.0.15
pipeline {
agent {label 'merlin'}
    stages {
        stage ("build-local-bwa-flow") {
            steps {
                 dir("ws-bwa-flow") {
		 checkout scm
                 script {
//                    git branch: 'release', url: 'git@github.com:falcon-computing/bwa-flow.git'
                        sh "rm -rf release"
                        sh "mkdir release"
                        dir("release"){
//                        sh "rsync -av --exclude=.* /curr/limark/genome-release/build/local/ /curr/limark/falcon2/"
//                        sh "rsync -av --exclude=.* /curr/limark/genome-release/build/common/ /curr/limark/falcon2/"
                        sh "source /curr/software/util/modules-tcl/init/bash"
                        sh "module load sdx/17.4; cmake -DCMAKE_BUILD_TYPE=Release -DRELEASE_VERSION=Internal on AWS -DDEPLOYMENT_DST=aws -DCMAKE_INSTALL_PREFIX=/curr/limark/falcon2/tools/bin .."
 //                       sh "cmake -DCMAKE_BUILD_TYPE=Release -DRELEASE_VERSION=aws -DDEPLOYMENT_DST= -DCMAKE_INSTALL_PREFIX=/curr/limark/falcon2/tools/bin .."
                        sh "make -j 8"
                        sh "make install"
			link = sh(returnStdout: true, script: 'cd /curr/limark/falcon2/tools/bin; link=s3://fcs-cicd-test/release/aws/bwa-flow/bwa-flow; echo $link; echo $link > latest')
                        sh "cd /curr/limark/falcon2/tools/bin; aws s3 cp bwa-flow s3://fcs-cicd-test/release/aws/bwa-flow/bwa-flow"
                        sh "cd /curr/limark/falcon2/tools/bin; aws s3 cp latest s3://fcs-cicd-test/release/aws/bwa-flow/latest"
                        sh "cd /curr/limark/falcon2/bin; rm -f latest"
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
	
