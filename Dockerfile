FROM ubuntu:18.04 AS build-env

WORKDIR /app


COPY . ./


RUN apt-get update && apt-get install -y liblapack-dev
RUN apt-get update && apt-get install -y cmake libboost-all-dev libarmadillo-dev libsnappy-dev
RUN rm -rf build ; mkdir build && cd build && cmake .. && make && make install

ENTRYPOINT ["ctmo", "test/Simulations/holstein_DMFT_T0.3/params1.json"]
