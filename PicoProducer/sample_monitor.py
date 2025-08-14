#!/usr/bin/env python3

import os
import sys
import json
import argparse
import subprocess
import re
from datetime import datetime
from collections import defaultdict

class SampleMonitor:
    """Monitor specific samples for TES, LTF, and JTF variations"""

    def __init__(self, year="2024", channel="mutau", verbose=False, config_file="sample_status.json"):
        self.year = year
        self.channel = channel
        self.verbose = verbose
        self.config_file = config_file
        
        # Define all specific samples based on your input
        self.samples = {
            "DY": [
                # DY to electrons
                "DYto2E_Bin-MLL-10to50",
                "DYto2E_Bin-MLL-50to120", 
                "DYto2E_Bin-MLL-120to200",
                "DYto2E_Bin-MLL-200to400",
                "DYto2E_Bin-MLL-400to800",
                "DYto2E_Bin-MLL-800to1500",
                "DYto2E_Bin-MLL-1500to2500",
                "DYto2E_Bin-MLL-2500to4000",
                "DYto2E_Bin-MLL-4000to6000",
                "DYto2E_Bin-MLL-6000",
                # DY to muons
                "DYto2Mu_Bin-MLL-10to50",
                "DYto2Mu_Bin-MLL-50to120",
                "DYto2Mu_Bin-MLL-120to200", 
                "DYto2Mu_Bin-MLL-200to400",
                "DYto2Mu_Bin-MLL-400to800",
                "DYto2Mu_Bin-MLL-800to1500",
                "DYto2Mu_Bin-MLL-1500to2500",
                "DYto2Mu_Bin-MLL-2500to4000",
                "DYto2Mu_Bin-MLL-4000to6000",
                "DYto2Mu_Bin-MLL-6000",
                # DY to taus
                "DYto2Tau_Bin-MLL-10to50",
                "DYto2Tau_Bin-MLL-50to120",
                "DYto2Tau_Bin-MLL-120to200",
                "DYto2Tau_Bin-MLL-200to400", 
                "DYto2Tau_Bin-MLL-400to800",
                "DYto2Tau_Bin-MLL-800to1500",
                "DYto2Tau_Bin-MLL-1500to2500",
                "DYto2Tau_Bin-MLL-2500to4000",
                "DYto2Tau_Bin-MLL-4000to6000",
                "DYto2Tau_Bin-MLL-6000"
            ],
            "TT": [
                "TTto2L2Nu",
                "TTtoLNu2Q", 
                "TTto4Q"
            ],
            "WJ": [
                "WtoTauNu-2Jets",
                "WtoMuNu-2Jets",
                "WtoENu-2Jets"
            ]
        }
        
        # Parse TES values from 0.97 to 1.03 with 0.002 step
        #self.tes_values = [round(0.97 + i * 0.002, 3) for i in range(31)]
        self.tes_values = [round(0.97 + i * 0.002, 3) for i in range(16)]  # 0.97 to 1.00 inclusive
        self.ltf_values = []
        self.jtf_values = []
        
        self.variations = {
            "tes": self.tes_values,
            "ltf": self.ltf_values, 
            "jtf": self.jtf_values
        }
    
    def get_all_samples(self, var_type):
        """Get sample list based on variation type"""
        if var_type == "jtf":
            # JTF uses DY + TT + WJ
            return self.samples["DY"] + self.samples["TT"] + self.samples["WJ"]
        else:
            # TES and LTF use DY + TT
            return self.samples["DY"] + self.samples["TT"]
    
    def create_tag(self, var_type, value):
        """Create tag string for pico.py commands"""
        value_str = f"{value:.3f}".replace('.', 'p')
        return f"_{var_type.upper()}{value_str}"
    
    def load_config(self):
        """Load status configuration from JSON"""
        if self.verbose:
            print(f"[VERBOSE] Loading config from: {self.config_file}")
        
        if os.path.exists(self.config_file):
            try:
                with open(self.config_file, 'r') as f:
                    config = json.load(f)
                if self.verbose:
                    print(f"[VERBOSE] Config loaded successfully")
                    print(f"[VERBOSE] Config contains {len(config.get('variations', {}))} variation types")
                    print(f"[VERBOSE] Last updated: {config.get('last_updated', 'Unknown')}")
                return config
            except Exception as e:
                if self.verbose:
                    print(f"[VERBOSE] Error loading config: {e}")
                    print(f"[VERBOSE] Initializing new config")
                return self.initialize_config()
        else:
            if self.verbose:
                print(f"[VERBOSE] Config file does not exist, initializing new config")
            return self.initialize_config()
    
    def save_config(self, config):
        """Save status configuration to JSON"""
        config["last_updated"] = datetime.now().isoformat()
        
        if self.verbose:
            print(f"[VERBOSE] Saving config to: {self.config_file}")
            total_jobs = sum(config.get('summary', {}).values())
            print(f"[VERBOSE] Config contains {total_jobs} total job entries")
        
        try:
            with open(self.config_file, 'w') as f:
                json.dump(config, f, indent=2)
            if self.verbose:
                print(f"[VERBOSE] Config saved successfully at {config['last_updated']}")
        except Exception as e:
            print(f"ERROR: Failed to save config: {e}")
            if self.verbose:
                print(f"[VERBOSE] Save error details: {type(e).__name__}: {str(e)}")
    
    def initialize_config(self):
        """Initialize configuration with all variations and samples"""
        config = {
            "variations": {},
            "last_updated": datetime.now().isoformat(),
            "summary": {
                "total_jobs": 0,
                "completed": 0,
                "failed": 0,
                "pending": 0
            }
        }
        
        total_jobs = 0
        
        for var_type, values in self.variations.items():
            config["variations"][var_type] = {}
            samples = self.get_all_samples(var_type)
            
            for value in values:
                value_key = f"{value:.3f}"
                config["variations"][var_type][value_key] = {}
                
                for sample in samples:
                    config["variations"][var_type][value_key][sample] = "pending"
                    total_jobs += 1
        
        config["summary"]["total_jobs"] = total_jobs
        self.save_config(config)
        return config
    
    def run_pico_command(self, action, var_type, value, sample):
        """Run pico.py command and return status"""
        tag = self.create_tag(var_type, value)
        
        if action == "status":
            command = f"pico.py status -y {self.year} -c {self.channel} -t {tag} -s {sample} -E {var_type}={value:.3f}"
        elif action == "submit":
            command = f"pico.py submit -y {self.year} -c {self.channel} -t {tag} -s {sample} -E {var_type}={value:.3f}"
        elif action == "resubmit":
            command = f"pico.py resubmit -y {self.year} -c {self.channel} -t {tag} -s {sample} -E {var_type}={value:.3f}"
        elif action == "hadd":
            # For hadd, we need all samples for this variation
            all_samples = self.get_all_samples(var_type)
            sample_str = " ".join(all_samples)
            command = f"pico.py hadd -y {self.year} -c {self.channel} -t {tag} -s {sample_str} -E {var_type}={value:.3f}"
        else:
            return None, "Unknown action"
        
        if self.verbose:
            print(f"[VERBOSE] Executing command: {command}")
            print(f"[VERBOSE] Working directory: {os.getcwd()}")
            print(f"[VERBOSE] Action: {action}, Sample: {sample}, Variation: {var_type}={value:.3f}")
        
        try:
            if self.verbose:
                print(f"[VERBOSE] Starting subprocess with timeout=36000s...")

            result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=36000)
            
            if self.verbose:
                print(f"[VERBOSE] Command completed with return code: {result.returncode}")
                print(f"[VERBOSE] STDOUT length: {len(result.stdout)} chars")
                print(f"[VERBOSE] STDERR length: {len(result.stderr)} chars")
                if result.stdout:
                    print(f"[VERBOSE] STDOUT: {result.stdout[:500]}{'...' if len(result.stdout) > 500 else ''}")
                if result.stderr:
                    print(f"[VERBOSE] STDERR: {result.stderr[:500]}{'...' if len(result.stderr) > 500 else ''}")
            
            return result, None
        except subprocess.TimeoutExpired:
            if self.verbose:
                print(f"[VERBOSE] Command timed out after 3600 seconds")
            return None, "Command timed out"
        except Exception as e:
            if self.verbose:
                print(f"[VERBOSE] Exception occurred: {type(e).__name__}: {str(e)}")
            return None, str(e)
    
    def parse_status(self, output):
        """Parse pico.py status output to determine job status"""
        if self.verbose:
            print(f"[VERBOSE] Parsing status from output (length: {len(output) if output else 0} chars)")
            if output:
                print(f"[VERBOSE] Raw output: '{output[:200]}{'...' if len(output) > 200 else ''}'")
        
        if not output or not output.strip():
            if self.verbose:
                print(f"[VERBOSE] Output is empty or whitespace only -> returning 'unknown'")
            return "unknown"
        
        output = output.upper()
        
        if "SUCCESS" in output and "FAIL" not in output and "MISS" not in output and "PEND" not in output:
            status = "completed"
        elif "FAIL" in output or "MISS" in output or "ERROR" in output:
            status = "failed"
        elif "PEND" in output or "SUBM" in output or "RUNNING" in output:
            status = "pending"
        else:
            status = "unknown"
        
        if self.verbose:
            print(f"[VERBOSE] Parsed status: '{status}' (based on keywords in output)")
        
        return status
    
    def check_status(self, var_types=None, auto_resubmit=False, auto_hadd=False):
        """Check status of all or specific variation types, and immediately resubmit/hadd and update config as soon as a single file is failed/completed."""
        if var_types is None:
            var_types = list(self.variations.keys())

        config = self.load_config()

        print("=== CHECKING STATUS ===")
        if auto_resubmit:
            print("(Auto-resubmit enabled - will resubmit failed jobs immediately)")
        if auto_hadd:
            print("(Auto-hadd enabled - will hadd completed variations immediately)")

        if self.verbose:
            print(f"[VERBOSE] Checking variation types: {var_types}")
            print(f"[VERBOSE] Config file: {self.config_file}")
            print(f"[VERBOSE] Year: {self.year}, Channel: {self.channel}")
            print(f"[VERBOSE] Auto-resubmit: {auto_resubmit}")
            print(f"[VERBOSE] Auto-hadd: {auto_hadd}")

        summary = {"completed": 0, "failed": 0, "pending": 0, "unknown": 0, "hadded": 0}
        total_resubmitted = 0
        total_checked = 0
        total_skipped = 0
        total_hadded = 0

        for var_type in var_types:
            if var_type not in self.variations:
                if self.verbose:
                    print(f"[VERBOSE] Skipping unknown variation type: {var_type}")
                continue

            print(f"\n--- {var_type.upper()} Variations ---")
            if self.verbose:
                print(f"[VERBOSE] Processing {len(self.variations[var_type])} values for {var_type}")

            for value in self.variations[var_type]:
                value_key = f"{value:.3f}"
                print(f"\n{var_type.upper()}={value:.3f}:")

                samples = self.get_all_samples(var_type)
                if self.verbose:
                    print(f"[VERBOSE] Checking {len(samples)} samples for {var_type}={value:.3f}")

                for i, sample in enumerate(samples):
                    if self.verbose:
                        print(f"[VERBOSE] [{i+1}/{len(samples)}] Checking sample: {sample}")

                    # Check if sample is already completed/hadded - skip status check if so
                    current_status = None
                    if (var_type in config["variations"] and 
                        value_key in config["variations"][var_type] and 
                        sample in config["variations"][var_type][value_key]):
                        current_status = config["variations"][var_type][value_key][sample]

                    # Skip checking if already in final state
                    if current_status in ["completed", "hadded"]:
                        status = current_status
                        print(f"  {sample}: {status.upper()} (cached)")
                        if self.verbose:
                            print(f"[VERBOSE] Skipping status check - already {status}")
                        total_skipped += 1
                    else:
                        # Perform actual status check
                        result, error = self.run_pico_command("status", var_type, value, sample)
                        total_checked += 1

                        if error:
                            status = "unknown"
                            print(f"  {sample}: ERROR - {error}")
                            if self.verbose:
                                print(f"[VERBOSE] Command failed for {sample}: {error}")
                        else:
                            status = self.parse_status(result.stdout)
                            print(f"  {sample}: {status.upper()}", end="")
                            if self.verbose:
                                print(f" [return_code={result.returncode}]", end="")

                        # Update config with current status immediately
                        if var_type in config["variations"] and value_key in config["variations"][var_type]:
                            old_status = config["variations"][var_type][value_key].get(sample, "unknown")
                            config["variations"][var_type][value_key][sample] = status
                            if self.verbose and old_status != status:
                                print(f"[VERBOSE] Status changed for {sample}: {old_status} → {status}")
                            self.save_config(config)

                        summary[status] += 1

                        # IMMEDIATE AUTO-RESUBMIT if failed and auto_resubmit is enabled
                        if status == "failed" and auto_resubmit:
                            print(" → Resubmitting...", end="", flush=True)
                            if self.verbose:
                                print(f"\n[VERBOSE] Auto-resubmitting failed job: {sample}")

                            resubmit_result, resubmit_error = self.run_pico_command("resubmit", var_type, value, sample)

                            if resubmit_error:
                                print(f" FAILED ({resubmit_error})")
                                status = "failed"  # Keep as failed
                                if self.verbose:
                                    print(f"[VERBOSE] Resubmit failed: {resubmit_error}")
                            elif resubmit_result.returncode == 0:
                                print(" SUCCESS")
                                status = "pending"  # Change to pending
                                total_resubmitted += 1
                                if self.verbose:
                                    print(f"[VERBOSE] Resubmit successful, status changed to pending")
                            else:
                                print(" FAILED (non-zero return code)")
                                status = "failed"  # Keep as failed
                                if self.verbose:
                                    print(f"[VERBOSE] Resubmit failed with return code: {resubmit_result.returncode}")
                            # Update config after resubmit
                            if var_type in config["variations"] and value_key in config["variations"][var_type]:
                                config["variations"][var_type][value_key][sample] = status
                                self.save_config(config)
                        elif status == "failed":
                            print()  # Just newline if not auto-resubmitting
                        else:
                            print()  # Newline for other statuses

                        # IMMEDIATE AUTO-HADD if completed and auto_hadd is enabled
                        if status == "completed" and auto_hadd:
                            print(" → HADDING...", end="", flush=True)
                            if self.verbose:
                                print(f"\n[VERBOSE] Auto-hadding for {var_type}={value:.3f} after {sample} completed")
                            # Only hadd for this variation (all samples)
                            result_hadd, error_hadd = self.run_pico_command("hadd", var_type, value, "")
                            if error_hadd:
                                print(f" FAILED ({error_hadd})")
                                if self.verbose:
                                    print(f"[VERBOSE] Hadd failed: {error_hadd}")
                            elif result_hadd.returncode == 0:
                                print(" SUCCESS")
                                total_hadded += 1
                                # Update all samples for this variation to "hadded"
                                all_samples = self.get_all_samples(var_type)
                                for s in all_samples:
                                    old_status = config["variations"][var_type][value_key].get(s, "unknown")
                                    config["variations"][var_type][value_key][s] = "hadded"
                                    if self.verbose:
                                        print(f"[VERBOSE] Status updated for {s}: {old_status} → hadded")
                                self.save_config(config)
                            else:
                                print(f" FAILED (return code: {result_hadd.returncode})")
                                if self.verbose:
                                    print(f"[VERBOSE] Hadd failed with return code: {result_hadd.returncode}")
                                    if result_hadd.stderr:
                                        print(f"[VERBOSE] Hadd stderr: {result_hadd.stderr[:200]}{'...' if len(result_hadd.stderr) > 200 else ''}")

        # Print summary
        print(f"\n=== SUMMARY ===")
        total_jobs = sum(summary.values())
        print(f"Total jobs: {total_jobs}")
        print(f"Completed: {summary['completed']}")
        print(f"Failed: {summary['failed']}")
        print(f"Pending: {summary['pending']}")
        print(f"Unknown: {summary['unknown']}")
        print(f"Hadded: {summary['hadded']}")
        print(f"Checked: {total_checked}, Cached: {total_skipped}")

        if self.verbose:
            completion_rate = (summary['completed'] + summary['hadded']) / total_jobs * 100 if total_jobs > 0 else 0
            failure_rate = summary['failed'] / total_jobs * 100 if total_jobs > 0 else 0
            cache_efficiency = total_skipped / (total_checked + total_skipped) * 100 if (total_checked + total_skipped) > 0 else 0
            print(f"[VERBOSE] Total completion rate (completed+hadded): {completion_rate:.1f}%")
            print(f"[VERBOSE] Failure rate: {failure_rate:.1f}%")
            print(f"[VERBOSE] Cache efficiency: {cache_efficiency:.1f}% (saved {total_skipped} status checks)")
            print(f"[VERBOSE] Config last updated: {config.get('last_updated', 'Unknown')}")

        return config
    
    def submit_variations(self, var_types=None, value_filter=None):
        """Submit jobs for variations"""
        if var_types is None:
            var_types = list(self.variations.keys())
        
        config = self.load_config()
        
        print("=== SUBMITTING VARIATIONS ===")
        if self.verbose:
            print(f"[VERBOSE] Submitting variation types: {var_types}")
            print(f"[VERBOSE] Value filter: {value_filter}")
        
        total_submitted = 0
        total_attempted = 0
        
        for var_type in var_types:
            if var_type not in self.variations:
                if self.verbose:
                    print(f"[VERBOSE] Skipping unknown variation type: {var_type}")
                continue
                
            print(f"\n--- Submitting {var_type.upper()} ---")
            if self.verbose:
                print(f"[VERBOSE] Processing {len(self.variations[var_type])} values for {var_type}")
            
            for value in self.variations[var_type]:
                if value_filter and value not in value_filter:
                    if self.verbose:
                        print(f"[VERBOSE] Skipping {var_type}={value:.3f} (not in filter)")
                    continue
                    
                print(f"\nSubmitting {var_type.upper()}={value:.3f}")
                samples = self.get_all_samples(var_type)
                if self.verbose:
                    print(f"[VERBOSE] Submitting {len(samples)} samples for {var_type}={value:.3f}")
                
                for i, sample in enumerate(samples):
                    if self.verbose:
                        print(f"[VERBOSE] [{i+1}/{len(samples)}] Submitting sample: {sample}")
                    
                    total_attempted += 1
                    result, error = self.run_pico_command("submit", var_type, value, sample)
                    
                    if error:
                        print(f"  {sample}: ERROR - {error}")
                        if self.verbose:
                            print(f"[VERBOSE] Submit failed for {sample}: {error}")
                    elif result.returncode == 0:
                        print(f"  {sample}: Submitted")
                        total_submitted += 1
                        # Update status in config
                        value_key = f"{value:.3f}"
                        if var_type in config["variations"] and value_key in config["variations"][var_type]:
                            old_status = config["variations"][var_type][value_key].get(sample, "unknown")
                            config["variations"][var_type][value_key][sample] = "pending"
                            if self.verbose:
                                print(f"[VERBOSE] Status updated for {sample}: {old_status} → pending")
                    else:
                        print(f"  {sample}: Failed to submit (return code: {result.returncode})")
                        if self.verbose:
                            print(f"[VERBOSE] Submit failed for {sample} with return code: {result.returncode}")
                            if result.stderr:
                                print(f"[VERBOSE] Submit stderr: {result.stderr[:200]}{'...' if len(result.stderr) > 200 else ''}")
        
        if self.verbose:
            print(f"[VERBOSE] Submission summary: {total_submitted}/{total_attempted} jobs submitted successfully")
            print(f"[VERBOSE] Saving config to {self.config_file}")
        
        self.save_config(config)
    
    def resubmit_failed(self, var_types=None):
        """Resubmit failed jobs"""
        if var_types is None:
            var_types = list(self.variations.keys())
        
        config = self.load_config()
        
        print("=== RESUBMITTING FAILED JOBS ===")
        
        resubmitted = 0
        
        for var_type in var_types:
            if var_type not in config["variations"]:
                continue
                
            print(f"\n--- Resubmitting {var_type.upper()} ---")
            
            for value_key, samples_dict in config["variations"][var_type].items():
                value = float(value_key)
                failed_samples = [sample for sample, status in samples_dict.items() if status == "failed"]
                
                if failed_samples:
                    print(f"\nResubmitting {var_type.upper()}={value:.3f} ({len(failed_samples)} failed samples)")
                    
                    for sample in failed_samples:
                        result, error = self.run_pico_command("resubmit", var_type, value, sample)
                        
                        if error:
                            print(f"  {sample}: ERROR - {error}")
                        elif result.returncode == 0:
                            print(f"  {sample}: Resubmitted")
                            config["variations"][var_type][value_key][sample] = "pending"
                            resubmitted += 1
                        else:
                            print(f"  {sample}: Failed to resubmit")
        
        print(f"\nTotal resubmitted: {resubmitted}")
        self.save_config(config)
    
    def hadd_completed(self, var_types=None):
        """Hadd completed variations"""
        if var_types is None:
            var_types = list(self.variations.keys())
        
        config = self.load_config()
        
        print("=== HADDING COMPLETED VARIATIONS ===")
        
        for var_type in var_types:
            if var_type not in config["variations"]:
                continue
                
            print(f"\n--- Hadding {var_type.upper()} ---")
            
            for value_key, samples_dict in config["variations"][var_type].items():
                value = float(value_key)
                
                # Check if all samples are completed
                all_samples = self.get_all_samples(var_type)
                completed_samples = [sample for sample in all_samples 
                                   if samples_dict.get(sample) == "completed"]
                
                if len(completed_samples) == len(all_samples):
                    print(f"\nHadding {var_type.upper()}={value:.3f} (all {len(all_samples)} samples completed)")
                    
                    result, error = self.run_pico_command("hadd", var_type, value, "")
                    
                    if error:
                        print(f"  ERROR: {error}")
                    elif result.returncode == 0:
                        print(f"  SUCCESS: Hadded {var_type.upper()}={value:.3f}")
                        # Mark as hadded
                        for sample in all_samples:
                            config["variations"][var_type][value_key][sample] = "hadded"
                    else:
                        print(f"  FAILED: Hadd failed for {var_type.upper()}={value:.3f}")
                elif completed_samples:
                    print(f"  {var_type.upper()}={value:.3f}: {len(completed_samples)}/{len(all_samples)} completed")
        
        self.save_config(config)
    
    def print_detailed_status(self):
        """Print detailed status from config file"""
        config = self.load_config()
        
        print("=== DETAILED STATUS REPORT ===")
        print(f"Last updated: {config.get('last_updated', 'Unknown')}")
        
        summary = config.get('summary', {})
        print(f"\nSUMMARY:")
        print(f"  Total jobs: {summary.get('total_jobs', 0)}")
        print(f"  Completed: {summary.get('completed', 0)}")
        print(f"  Failed: {summary.get('failed', 0)}")
        print(f"  Pending: {summary.get('pending', 0)}")
        print(f"  Hadded: {summary.get('hadded', 0)}")
        
        for var_type, variations in config.get('variations', {}).items():
            print(f"\n--- {var_type.upper()} DETAILS ---")
            
            for value_key, samples_dict in variations.items():
                statuses = list(samples_dict.values())
                completed = statuses.count('completed')
                failed = statuses.count('failed') 
                pending = statuses.count('pending')
                hadded = statuses.count('hadded')
                total = len(statuses)
                
                print(f"  {var_type.upper()}={value_key}: {completed}C/{failed}F/{pending}P/{hadded}H (total: {total})")


def main():
    parser = argparse.ArgumentParser(description="Sample status monitor for TES/LTF/JTF variations")
    
    # Basic options
    parser.add_argument("-y", "--year", default="2024", help="Year (default: 2024)")
    parser.add_argument("-c", "--channel", default="mutau", help="Channel (default: mutau)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--config", default="sample_status.json", help="JSON config file (default: sample_status.json)")
    parser.add_argument("--auto-resubmit", action="store_true", help="Automatically resubmit failed jobs during status check (default: False)")
    parser.add_argument("--auto-hadd", action="store_true", help="Automatically hadd completed variations during status check (default: False)")
    
    # Actions
    parser.add_argument("action", choices=["status", "submit", "resubmit", "hadd", "report"], 
                       help="Action to perform")
    
    # Filters
    parser.add_argument("--var-types", nargs="*", choices=["tes", "ltf", "jtf"], 
                       help="Variation types to process")
    parser.add_argument("--values", nargs="*", type=float,
                       help="Specific variation values to process")
    
    args = parser.parse_args()
    
    # Initialize monitor
    monitor = SampleMonitor(
        year=args.year,
        channel=args.channel,
        verbose=args.verbose,
        config_file=args.config
    )
    
    if args.verbose:
        print(f"[VERBOSE] ======= SAMPLE MONITOR STARTED =======")
        print(f"[VERBOSE] Action: {args.action}")
        print(f"[VERBOSE] Year: {args.year}")
        print(f"[VERBOSE] Channel: {args.channel}")
        print(f"[VERBOSE] Config file: {args.config}")
        print(f"[VERBOSE] Auto-resubmit: {args.auto_resubmit}")
        print(f"[VERBOSE] Auto-hadd: {args.auto_hadd}")
        print(f"[VERBOSE] Variation types filter: {args.var_types}")
        print(f"[VERBOSE] Values filter: {args.values}")
        print(f"[VERBOSE] Working directory: {os.getcwd()}")
        print(f"[VERBOSE] Python executable: {sys.executable}")
        print(f"[VERBOSE] Script path: {__file__}")
        
        # Print variation summary
        for var_type, values in monitor.variations.items():
            samples = monitor.get_all_samples(var_type)
            total_jobs = len(values) * len(samples)
            print(f"[VERBOSE] {var_type.upper()}: {len(values)} values × {len(samples)} samples = {total_jobs} jobs")
        
        print(f"[VERBOSE] ========================================")
    
    # Execute action
    try:
        if args.action == "status":
            monitor.check_status(args.var_types, auto_resubmit=args.auto_resubmit, auto_hadd=args.auto_hadd)
        elif args.action == "submit":
            monitor.submit_variations(args.var_types, args.values)
        elif args.action == "resubmit":
            monitor.resubmit_failed(args.var_types)
        elif args.action == "hadd":
            monitor.hadd_completed(args.var_types)
        elif args.action == "report":
            monitor.print_detailed_status()
        
        if args.verbose:
            print(f"[VERBOSE] Action '{args.action}' completed successfully")
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {str(e)}")
        if args.verbose:
            import traceback
            print(f"[VERBOSE] Full traceback:")
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
