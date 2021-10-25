FROM osrf/ros:noetic-desktop-full-focal

ARG PKGS="sudo \
  ssh \
  openssh-client \
  git \
  build-essential \
  libgtest-dev \
  libglfw3-dev \
  libglu1-mesa-dev \
  vim \
  unzip \
  python3-catkin-tools \
  libsuitesparse-dev \
  freeglut3-dev \
  python3-pip"

RUN apt-get update && \
    apt-get dist-upgrade -y && \
    apt-get install -y ${PKGS}

RUN mkdir -p ${HOME}/workspace/src
RUN git clone https://github.com/srrg-sapienza/srrg2_solver.git ${HOME}/workspace/src/srrg2_solver
RUN git clone https://github.com/srrg-sapienza/srrg2_core.git ${HOME}/workspace/src/srrg2_core
RUN git clone https://github.com/srrg-sapienza/srrg2_cmake_modules.git ${HOME}/workspace/src/srrg_cmake_modules
RUN git clone https://github.com/srrg-sapienza/srrg2-hipe.git ${HOME}/workspace/src/srrg2_hipe

RUN python3 -m pip install osrf-pycommon
RUN cd ${HOME}/workspace/src/srrg2_solver && rm -rf srrg2_solver_gui srrg2_solver_extras srrg2_solver_experiments
RUN cd ${HOME}/workspace/src/srrg2_core && rm -rf srrg2_core_ros 

RUN cd ${HOME}/workspace/ && . /opt/ros/noetic/setup.sh && cd src && catkin_init_workspace && cd .. && catkin_make -j4 
RUN echo "source ${HOME}/workspace/devel/setup.bash" >> ${HOME}/.bashrc
WORKDIR /root/
