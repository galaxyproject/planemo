pipeline {
  agent {
    dockerfile {
      filename 'tests/Dockerfile.jenkins'
      additionalBuildArgs '--build-arg SKIP_TOX=true'
    }
  }

  stages {
    stage("test") {
      steps {
        sh 'tox -e py27-quick'
      }
    }
  }

}
