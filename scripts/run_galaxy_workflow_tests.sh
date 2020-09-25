#!/bin/bash

# Usage: http://planemo.readthedocs.io/en/latest/test_format.html#galaxy-testing-template

: ${PLANEMO_TARGET:="planemo==0.52.0"}
: ${PLANEMO_OPTIONS:=""}  # e.g. PLANEMO_OPTIONS="--verbose"
: ${PLANEMO_PROFILE_NAME:="wxflowtest"}
: ${PLANEMO_SERVE_PORT:="9019"}
: ${PLANEMO_GALAXY_BRANCH:="master"}
: ${PLANEMO_TEST_STYLE:="serve_and_test"}  # profile_serve_and_test, serve_and_test, docker_serve_and_test, test, docker_test, docker_test_path_paste
: ${PLANEMO_SERVE_DATABASE_TYPE:="postgres"}  # used if not using Docker with PLANEMO_TEST_STYLE
: ${PLANEMO_DOCKER_GALAXY_IMAGE:="quay.io/bgruening/galaxy:20.05"}  # used if used Docker with PLANEMO_TEST_STYLE
: ${PLANEMO_VIRTUAL_ENV:=".venv"}
: ${GALAXY_URL:="http://localhost:$PLANEMO_SERVE_PORT"}

# Ensure Planemo is installed.
if [ ! -d "${PLANEMO_VIRTUAL_ENV}" ]; then
    virtualenv "${PLANEMO_VIRTUAL_ENV}"
    . "${PLANEMO_VIRTUAL_ENV}"/bin/activate
    pip install -U pip>7
    # Intentionally expand wildcards in PLANEMO_TARGET.
    shopt -s extglob
    pip install ${PLANEMO_TARGET}
fi
. "${PLANEMO_VIRTUAL_ENV}"/bin/activate

# Run test.
# This example shows off a bunch of different ways one could test with Planemo,
# but for actual workflow testing projects - probably best just to take one of the last
# two very easy invocations to simplify things.
if [ "$PLANEMO_TEST_STYLE" = "profile_serve_and_test" ]; then
    planemo $PLANEMO_OPTIONS profile_create \
        --database_type "$PLANEMO_SERVE_DATABASE_TYPE" \
        "$PLANEMO_PROFILE_NAME"
    planemo $PLANEMO_OPTIONS serve \
        --daemon \
        --galaxy_branch "$PLANEMO_GALAXY_BRANCH" \
        --profile "$PLANEMO_PROFILE_NAME" \
        --port "$PLANEMO_SERVE_PORT" \
        "$1"
    planemo $PLANEMO_OPTIONS test \
        --galaxy_url "$GALAXY_URL" \
        --engine external_galaxy \
        "$1"
elif [ "$PLANEMO_TEST_STYLE" = "serve_and_test" ]; then
    planemo $PLANEMO_OPTIONS serve \
        --daemon \
        --galaxy_branch "$PLANEMO_GALAXY_BRANCH" \
        --database_type "$PLANEMO_SERVE_DATABASE_TYPE" \
        --port "$PLANEMO_SERVE_PORT" \
        "$1"
    planemo $PLANEMO_OPTIONS test \
        --galaxy_url "$GALAXY_URL" \
        --engine external_galaxy \
        "$1"
elif [ "$PLANEMO_TEST_STYLE" = "docker_serve_and_test" ]; then
    docker pull "${PLANEMO_DOCKER_GALAXY_IMAGE}"
    planemo $PLANEMO_OPTIONS serve \
        --daemon \
        --engine docker_galaxy \
        --docker_galaxy_image "${PLANEMO_DOCKER_GALAXY_IMAGE}" \
        --port "$PLANEMO_SERVE_PORT" \
        "$1"
    planemo $PLANEMO_OPTIONS test \
        --galaxy_url "$GALAXY_URL" \
        --engine external_galaxy \
        "$1"
elif [ "$PLANEMO_TEST_STYLE" = "test" ]; then
    # TODO: this conda_init shouldn't be needed, but this mode is broken without it.
    planemo conda_init || true

    planemo $PLANEMO_OPTIONS test \
        --database_type "$PLANEMO_SERVE_DATABASE_TYPE" \
        --galaxy_branch "$PLANEMO_GALAXY_BRANCH" \
        "$1"
elif [ "$PLANEMO_TEST_STYLE" = "docker_test" ]; then
    # TODO: This variant isn't super usable yet because there is too much logging, hence the dev null
    # redirect.
    docker pull "${PLANEMO_DOCKER_GALAXY_IMAGE}"
    planemo $PLANEMO_OPTIONS test \
        --engine docker_galaxy \
        --docker_galaxy_image "${PLANEMO_DOCKER_GALAXY_IMAGE}" \
        "$1" > /dev/null
elif [ "$PLANEMO_TEST_STYLE" = "docker_test_path_paste" ]; then
    # Same as above but mount the test data and use file:// path pastes when uploading
    # files (more robust and quick if working with really large files).
    docker pull "${PLANEMO_DOCKER_GALAXY_IMAGE}"
    planemo $PLANEMO_OPTIONS test \
        --engine docker_galaxy \
        --docker_extra_volume . \
        --paste_test_data_paths \
        --docker_galaxy_image "${PLANEMO_DOCKER_GALAXY_IMAGE}" \
        "$1" > /dev/null
elif [ "$PLANEMO_TEST_STYLE" = "manual_docker_run_and_test" ]; then
    docker pull "${PLANEMO_DOCKER_GALAXY_IMAGE}"
    docker run -d -e "NONUSE=nodejs,proftp,reports" -p "${PLANEMO_SERVE_PORT}:80" "${PLANEMO_DOCKER_GALAXY_IMAGE}"
    galaxy-wait -g "http://localhost:${PLANEMO_SERVE_PORT}"
    planemo $PLANEMO_OPTIONS test \
        --engine external_galaxy \
        --galaxy_url "$GALAXY_URL" \
        --galaxy_admin_key admin \
        --galaxy_user_key admin \
        "$1"
elif [ "$PLANEMO_TEST_STYLE" = "external_galaxy" ]; then
    if [[ -n $PLANEMO_INSTALL_TOOLS ]]; then
        INSTALL_TOOLS="";
    else
        INSTALL_TOOLS="--no_shed_install";
    fi
    planemo $PLANEMO_OPTIONS test \
        --engine external_galaxy \
        --galaxy_url "$GALAXY_URL" \
        --galaxy_admin_key "$PLANEMO_ADMIN_KEY" \
        --galaxy_user_key "$PLANEMO_USER_KEY" \
        $INSTALL_TOOLS \
        "$1"
else
    echo "Unknown test style ${PLANEMO_TEST_STYLE}"
    exit 1
fi
