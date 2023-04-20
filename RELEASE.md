# How to release this project?

This project follows the gitflow model, so we have `develop` as development branch and `master` as production branch.

In order to release this project, you must execute from `develop` (with latest changes) the following command:
```commandline
mvn -B gitflow:release-start
```
Then, all the suffixes `-SNAPSHOT` will be removed and the `release/X.X.X` branch will be created, so then you must execute:
```commandline
mvn -B gitflow:release-finish
```
So in develop, the version of the poms will be incremented automatically and master will contain the release branch merged and the tag created.

For more information about how it works the plugin (options, parameters, etc) please refer to this [link](https://github.com/aleksandr-m/gitflow-maven-plugin).