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
        self.tes_values = [round(0.97 + i * 0.002, 3) for i in range(31)]
        self.ltf_values = [0.97, 1.03]
        self.jtf_values = [0.90, 1.10]
        
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
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as f:
                return json.load(f)
        else:
            return self.initialize_config()
    
    def save_config(self, config):
        """Save status configuration to JSON"""
        config["last_updated"] = datetime.now().isoformat()
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=2)
    
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
            print(f"Running: {command}")
        
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=60)
            return result, None
        except subprocess.TimeoutExpired:
            return None, "Command timed out"
        except Exception as e:
            return None, str(e)
    
    def parse_status(self, output):
        """Parse pico.py status output to determine job status"""
        if not output or not output.strip():
            return "unknown"
        
        output = output.upper()
        
        if "SUCCESS" in output and "FAIL" not in output and "MISS" not in output and "PEND" not in output:
            return "completed"
        elif "FAIL" in output or "MISS" in output or "ERROR" in output:
            return "failed"
        elif "PEND" in output or "SUBM" in output or "RUNNING" in output:
            return "pending"
        else:
            return "unknown"
    
    def check_status(self, var_types=None, auto_resubmit=False):
        """Check status of all or specific variation types"""
        if var_types is None:
            var_types = list(self.variations.keys())
        
        config = self.load_config()
        
        print("=== CHECKING STATUS ===")
        if auto_resubmit:
            print("(Auto-resubmit enabled)")
        
        summary = {"completed": 0, "failed": 0, "pending": 0, "unknown": 0}
        failed_jobs = []  # Track failed jobs for auto-resubmit
        
        for var_type in var_types:
            if var_type not in self.variations:
                continue
                
            print(f"\n--- {var_type.upper()} Variations ---")
            
            for value in self.variations[var_type]:
                value_key = f"{value:.3f}"
                print(f"\n{var_type.upper()}={value:.3f}:")
                
                samples = self.get_all_samples(var_type)
                
                for sample in samples:
                    result, error = self.run_pico_command("status", var_type, value, sample)
                    
                    if error:
                        status = "unknown"
                        print(f"  {sample}: ERROR - {error}")
                    else:
                        status = self.parse_status(result.stdout)
                        print(f"  {sample}: {status.upper()}")
                    
                    # Track failed jobs for auto-resubmit
                    if status == "failed" and auto_resubmit:
                        failed_jobs.append((var_type, value, sample))
                    
                    # Update config
                    if var_type in config["variations"] and value_key in config["variations"][var_type]:
                        config["variations"][var_type][value_key][sample] = status
                    
                    summary[status] += 1
        
        # Update summary in config
        config["summary"].update(summary)
        config["summary"]["total_jobs"] = sum(summary.values())
        
        self.save_config(config)
        
        # Auto-resubmit failed jobs if enabled
        if auto_resubmit and failed_jobs:
            print(f"\n=== AUTO-RESUBMITTING {len(failed_jobs)} FAILED JOBS ===")
            resubmitted = 0
            
            for var_type, value, sample in failed_jobs:
                result, error = self.run_pico_command("resubmit", var_type, value, sample)
                
                if error:
                    print(f"  {sample} ({var_type.upper()}={value:.3f}): ERROR - {error}")
                elif result.returncode == 0:
                    print(f"  {sample} ({var_type.upper()}={value:.3f}): Resubmitted")
                    # Update status in config
                    value_key = f"{value:.3f}"
                    config["variations"][var_type][value_key][sample] = "pending"
                    resubmitted += 1
                    # Update summary
                    summary["failed"] -= 1
                    summary["pending"] += 1
                else:
                    print(f"  {sample} ({var_type.upper()}={value:.3f}): Failed to resubmit")
            
            print(f"Auto-resubmitted: {resubmitted}/{len(failed_jobs)} jobs")
            
            # Update config with new statuses
            config["summary"].update(summary)
            self.save_config(config)
        
        # Print summary
        print(f"\n=== SUMMARY ===")
        print(f"Total jobs: {summary['completed'] + summary['failed'] + summary['pending'] + summary['unknown']}")
        print(f"Completed: {summary['completed']}")
        print(f"Failed: {summary['failed']}")
        print(f"Pending: {summary['pending']}")
        print(f"Unknown: {summary['unknown']}")
        
        return config
    
    def submit_variations(self, var_types=None, value_filter=None):
        """Submit jobs for variations"""
        if var_types is None:
            var_types = list(self.variations.keys())
        
        config = self.load_config()
        
        print("=== SUBMITTING VARIATIONS ===")
        
        for var_type in var_types:
            if var_type not in self.variations:
                continue
                
            print(f"\n--- Submitting {var_type.upper()} ---")
            
            for value in self.variations[var_type]:
                if value_filter and value not in value_filter:
                    continue
                    
                print(f"\nSubmitting {var_type.upper()}={value:.3f}")
                samples = self.get_all_samples(var_type)
                
                for sample in samples:
                    result, error = self.run_pico_command("submit", var_type, value, sample)
                    
                    if error:
                        print(f"  {sample}: ERROR - {error}")
                    elif result.returncode == 0:
                        print(f"  {sample}: Submitted")
                        # Update status in config
                        value_key = f"{value:.3f}"
                        if var_type in config["variations"] and value_key in config["variations"][var_type]:
                            config["variations"][var_type][value_key][sample] = "pending"
                    else:
                        print(f"  {sample}: Failed to submit")
        
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
    
    # Execute action
    if args.action == "status":
        monitor.check_status(args.var_types, auto_resubmit=args.auto_resubmit)
    elif args.action == "submit":
        monitor.submit_variations(args.var_types, args.values)
    elif args.action == "resubmit":
        monitor.resubmit_failed(args.var_types)
    elif args.action == "hadd":
        monitor.hadd_completed(args.var_types)
    elif args.action == "report":
        monitor.print_detailed_status()


if __name__ == "__main__":
    main()
