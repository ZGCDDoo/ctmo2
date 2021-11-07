FROM ubuntu:20.04 AS build-env

WORKDIR /app


COPY . ./

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential liblapack-dev cmake libboost-all-dev libarmadillo-dev libsnappy-dev
RUN rm -rf build ; mkdir build && cd build && cmake .. && make && make install

ENTRYPOINT ["ctmo"]
