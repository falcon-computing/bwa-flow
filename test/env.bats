#!/usr/bin/env bats
load global

@test "check bwa binary" {
  [ ! -z "$bwa_bin" ]
  [ -f "$bwa_bin" ]
}

@test "check sambamba" {
  [ ! -z "$sambamba" ]
  [ -f "$sambamba" ]
}

@test "check if temp dir is writable" {
  run mkdir -p $temp_dir
  [ -d $temp_dir ]
}

@test "check reference genome" {
  [ ! -z "$ref_genome" ]
  [ -f "$ref_genome" ]
  [ -f "${ref_genome}.pac" ]
  [ -f "${ref_genome}.sa" ]
}

