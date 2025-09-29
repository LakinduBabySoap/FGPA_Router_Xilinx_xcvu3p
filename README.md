# FPGA Router for Xilinx xcvu3p  

This repository contains a high-performance FPGA router targeting the Xilinx xcvu3p device. It was developed for the CENG 4120 course and solves the FPGA routing problem by generating congestion-free routing trees for given nets using a simplified routing resource graph.

---

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Problem Description](#problem-description)
- [Input & Output Formats](#input--output-formats)
- [Features](#features)
- [Algorithm Details](#algorithm-details)
- [Optimizations](#optimizations)
- [Installation](#installation)
- [Usage](#usage)
- [Benchmarks & Grading](#benchmarks--grading)
- [Limitations](#limitations)
- [FAQ](#faq)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Overview

FPGA routing is a critical phase in the FPGA compilation flow, where logic elements are interconnected using the device's routing resources. This router takes a netlist and a device graph as input and outputs a routing solution in the required format.

**Highlights:**
- **Target Device:** Xilinx xcvu3p
- **Language:** C++ (with OpenMP for parallelism)
- **Algorithm:** Bidirectional A* search for pathfinding
- **Optimizations:** Selective node loading, memory mapping, early termination, beam search
- **Performance Goal:** Route designs within strict time limits (100–250s on HPC servers)

---

## Background

FPGAs consist of logic resources (LUTs, FFs) and routing resources (wires, programmable interconnects). After placement, the router connects these resources based on the netlist, ensuring no shared nodes between nets (to avoid congestion).

**Key Concepts:**
- **Tile:** Basic grid unit (e.g., CLE, INT)
- **Node:** Electrically connected wire, directional, with coordinates
- **PIP:** Programmable Interconnect Point (directed edge)
- **Net:** Source node and multiple sinks to connect via a tree

---

## Problem Description

Given:
- **G = (V, E):** Routing resource graph (nodes and PIPs)
- **N:** Netlist (each net has a source and sinks)

**Goal:**  
For each net, generate a tree connecting source to sinks, with no node sharing across nets (no congestion).

---

## Input & Output Formats

### Input

1. **Device File** (`xcvu3p.device`):  
   - First line: Number of nodes `n`
   - Next `n` lines:  
     `<node ID> <node type> <node length> <begin x> <begin y> <end x> <end y> <node name>`
   - Remaining lines:  
     `<Parent Node ID> <Child Node ID 0> [<Child Node ID 1> ...]`

2. **Netlist File** (`<design>.netlist`):  
   - First line: Number of nets `m`
   - Next `m` lines:  
     `<net ID> <net name> <source node ID> [<sink node ID 0> ...]`

### Output

- **Route File** (`<design>.route`):
  ```
  <net ID> <net name>
  <Parent node1 ID> <Child node1 ID>
  <Parent node2 ID> <Child node2 ID>
  ...

  <net ID> <net name>
  ...
  ```

---

## Features

- **Efficient Parsing:** Single-pass parsing of large device files using memory mapping and parallel processing
- **Bounding Box Optimization:** Loads only nodes within an expanded bounding box around source/sink coordinates
- **Bidirectional A* Search:** Finds paths from source to each sink, avoiding used nodes
- **Congestion Handling:** Tracks used nodes to avoid congestion (no rip-up/re-route)
- **Parallelism:** Uses OpenMP for parsing and edge processing
- **Heuristic:** Manhattan distance for A* f-score
- **Early Termination & Beam Search:** Prunes search space for speed

---

## Algorithm Details

1. **Parsing:**
   - Load netlist into memory
   - Parse device file:  
     - Collect required nodes (sources/sinks)
     - Compute bounding box
     - Load nodes in box, then parse edges using mmap and OpenMP

2. **Routing:**
   - Track global `used_nodes` (sources/sinks initially)
   - For each net:
     - Route to each sink using bidirectional A*
     - Allow intra-net node reuse, avoid global used nodes
     - Add intermediate nodes to global used set after routing

3. **Bidirectional A* Search:**
   - Forward search from source, backward from sink
   - Priority queues with f-score = g + heuristic (Manhattan)
   - Meet at a point, reconstruct path as PIP pairs
   - Skip used nodes except target
   - Limits: Max iterations (10k), beam width (500), early stop after 800 stagnant iterations

4. **Output:**  
   - Write PIPs per net in the required format

---

## Optimizations

- **Memory Efficiency:** Skips "NULL" and "CONFIG" nodes; uses hash maps/sets for O(1) lookups
- **Performance:**  
  - Parallel node parsing (`#pragma omp parallel for`)
  - Memory-mapped edge parsing with chunked parallel processing
  - Bounding box expansion (1.3x) to focus on relevant graph subset
- **Search Pruning:** Early termination, beam width to handle large graphs
- **Thread Limiting:** Caps OpenMP threads to 8 for HPC compliance
- **Timing:** Measures total execution time

---

## Installation

### Dependencies

- C++11 or later compiler (e.g., g++)
- OpenMP (enabled by default in modern g++)
- Unix-like system (for mmap, sysconf; tested on Linux/HPC)
- To download the device file copy and paste the following link to your browser (if prompted for download, just click download since sometimes google asks for confirmation due to large file size)
```sh
https://drive.usercontent.google.com/download?id=1IHtqsgXixh6pMkBdi1LM3KuQWYV4fHmQ&export=download&authuser=0
```
-To extract the downloaded file, then open the terminal and type the following
```sh
tar -xjf xcvu3p.tar.bz2
```
- The xcvu3p.device file will then be created and for the router remember to enter the full location of the device file (type in pwd to get the current directory path)

### Build

```sh
g++ -std=c++11 -fopenmp -O3 project_finalver.cpp -o router
```

---

## Usage

```sh
./router <netlist_file> <output_file> <device_file>
```
*Note- Make sure to enter the full file paths to be safe (used pwd to get the directory path)

**Example:**
```sh
./router benchmarks/design3.netlist output/design3.route xcvu3p.device
```

- Loads nets from netlist
- Parses device selectively
- Routes using bidirectional A*
- Writes route file
- Prints progress and timing

**On HPC (the server used for testing project):**
- Use Slurm:  
  `srun -p hpc_72h -w hpc[11-14] --cpus-per-task=8 --pty bash -i`
- Test runtime to meet limits
Note - Number of threads is limited to 8 for the project!
---

## Benchmarks & Grading

- **Benchmarks:** 5 designs (1-2 easy, 3-4 medium, 5 real FPGA signals)
- **Time limits:** 100s (1-4), 250s (5) on HPC11-14
- **Ranking:** No congestion > Max routed nets > Min wirelength
- **Weighted sum:** (0.1, 0.1, 0.2, 0.2, 0.4)
- **Scores:** 100% (1st) to 50% (others)
- **Evaluator:**  
  `./eval <device> <netlist> <result>`

---

## Limitations

- **No Congestion Resolution:** Avoids used nodes but doesn't rip-up/re-route if paths fail
- **Single-Path per Sink:** Routes source to each sink independently; no shared tree branches within net (suboptimal wirelength)
- **Memory Intensive:** Still loads significant graph portion; may fail on very large designs without further pruning
- **No Multi-Sink Optimization:** Treats multi-sink nets as separate paths
- **Fallback on mmap Failure:** Falls back to traditional reading (not fully implemented)
- **Assumes Valid Input:** No extensive error checking

---

## FAQ

- **How to reduce device file reading time?**  
  Use multi-threading (implemented)

- **How to view large files?**  
  Use `less` or `tail -n <lines> <file>`.

---

## Acknowledgments

- Course instructor and TA (for the specifications)
- Xilinx/Vivado concepts
- OpenMP and C++ standard library

---

## Libraries Used

```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <climits>
#include <chrono>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
```