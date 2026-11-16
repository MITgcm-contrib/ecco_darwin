# Working on Pleiades with Claude Code

*NASA NAS · Pleiades · Claude Code*

A step-by-step guide for running commands, checking jobs, and copying files on NASA's Pleiades supercomputer by asking Claude Code on your Mac. Setup takes about 15 minutes. After that, each working day starts with one login.

| You need | Notes |
|---|---|
| A NAS account and your NAS username | The username you use for ssh to NAS (it may differ from your Mac username). |
| Your NAS two-factor token (RSA SecurID passcode + PIN) | The same login you already use to reach Pleiades. |
| A Mac with the Terminal app | Applications → Utilities → Terminal. |
| Claude Code installed | Check with `claude --version` in Terminal. If it is missing, install it from Anthropic's Claude Code documentation. |

Throughout this guide, replace `<username>` with your NAS username. Text in code blocks is typed exactly as shown (one line at a time), or is what you should see.

## Contents

- [How it works](#how-it-works)
- [One-time setup (about 15 minutes)](#one-time-setup-about-15-minutes)
- [Each working session](#each-working-session)
- [Waiting for long jobs](#waiting-for-long-jobs)
- [Security and data](#security-and-data)
- [Troubleshooting](#troubleshooting)
- [Quick reference](#quick-reference)

## How it works

NAS asks for your passcode and PIN at its secure front end, `sfe`, and every connection to Pleiades (`pfe`) passes through it. Normally that means typing your passcode for every ssh command, which Claude cannot do.

SSH **connection sharing** solves this. You log in once; SSH keeps that authenticated connection open in the background and saves a small socket file in `~/.ssh/`. Every later `ssh pfe` or `scp` command, whether you or Claude runs it, reuses the open connection and needs no passcode.

| Where | What runs there |
|---|---|
| Your Mac | Claude Code. It decides what to run, sends each `ssh` or `scp` command, and reads the output. |
| `pfe` (Pleiades front end) | Only the individual commands Claude sends (`qstat`, `ls`, `tail`, a short script). Each one exits when it finishes. Nothing is installed on Pleiades and no Claude program runs there. |
| Compute nodes | Only the PBS jobs you choose to submit. |
| Anthropic | The Claude model. It sees the output of each command as part of your conversation. |

To NAS this looks exactly like you logging in and typing the commands yourself: everything runs under your account and appears in the usual NAS logs.

## One-time setup (about 15 minutes)

### 1. Check that your Mac already reaches Pleiades

Open Terminal and type:

```bash
ssh pfe
```

Enter your passcode and PIN when asked. You should end up at a Pleiades prompt that shows a `pfe` host name. Type `exit` to come back to your Mac.

> **Note:** If `ssh pfe` gives `Could not resolve hostname pfe`, your Mac is not set up for NAS yet. Follow the NAS knowledge-base article "Setting Up the SSH Passthrough Feature" (search for it on the NAS HECC support site), then come back to this step. It creates `~/.ssh/config` with the `sfe` and `pfe` entries this guide builds on.

### 2. Make a backup of your SSH configuration file

```bash
cp ~/.ssh/config ~/.ssh/config.backup
```

Nothing is printed. If anything goes wrong later, `cp ~/.ssh/config.backup ~/.ssh/config` puts the original back.

### 3. Add the connection-sharing settings to the top of the file

Open the file in the nano text editor:

```bash
nano ~/.ssh/config
```

The cursor starts at the very top of the file. Press **Return** twice to make two empty lines, press the **up arrow** twice to go back to the top, and paste these four lines:

```sshconfig
Host sfe pfe lfe
  ControlMaster auto
  ControlPath ~/.ssh/cm-%r@%h:%p
  ControlPersist 12h
```

Save and close: press **Control-O**, then **Return** (save), then **Control-X** (exit).

> **Important:** These lines must be at the top. SSH uses the first setting it finds for each host, and the NAS entries further down the file would otherwise take precedence. Do not change the NAS entries below them.

What the lines do: `ControlMaster auto` shares one connection; `ControlPath` is where the socket file is kept; `ControlPersist 12h` keeps the connection open for 12 hours after the last command that used it.

### 4. Make the SSH folder private

```bash
chmod 700 ~/.ssh
chmod 600 ~/.ssh/config
```

Nothing is printed. This stops other users of your Mac from reading your SSH files or using your open connection.

### 5. Test the shared connection in Terminal

Open it (enter your passcode and PIN):

```bash
ssh -fN pfe
```

After the login, the command returns straight to your Mac prompt; the connection stays open in the background (`-f`) without opening a remote shell (`-N`). Now check it:

```bash
ssh -O check pfe
```

You should see:

```text
Master running (pid=12345)
```

Finally run a command on Pleiades. It should answer immediately, with no passcode prompt:

```bash
ssh -o BatchMode=yes pfe hostname
```

You should see:

```text
pfe23
```

Close the test connection with `ssh -O exit pfe`. Setup is done.

### 6. Optional: give Claude the NAS ground rules

Claude Code reads a file called `CLAUDE.md` at the start of every session. Put the rules below in `~/.claude/CLAUDE.md` (applies to all your projects) or in a `CLAUDE.md` file in your project folder. Create or open the file with `nano ~/.claude/CLAUDE.md`, paste, then **Control-O**, **Return**, **Control-X**.

````markdown
## Working on NASA NAS (Pleiades)
- Reach Pleiades only with `ssh -o BatchMode=yes pfe '<command>'` and
  `scp -o BatchMode=yes`. If that fails, ask me to run `! ssh -fN pfe`;
  never try to log in yourself or handle my passcode.
- Use full paths for PBS commands: /PBS/bin/qstat, /PBS/bin/qsub.
- Ask me before every qsub or qdel (jobs use my allocation).
- Treat my existing run directories as read-only; work in copies.
- Never kill processes or jobs by name pattern; only by ID, and only ones
  you started.
- Run analysis on pfe next to the data; copy back only small results.
- Never work around NAS security settings.
````

## Each working session

### 1. Start Claude Code in your project folder

```bash
cd ~/path/to/your/project
claude
```

### 2. Open the NAS connection from inside Claude Code

At the Claude Code prompt, type the line below. The `!` runs it as a normal shell command in your own terminal, so you type the passcode and PIN yourself; Claude never sees them.

```text
! ssh -fN pfe
```

When you are back at the Claude prompt, the connection is open.

### 3. Ask Claude for what you need, in plain English

For example:

- "Check my jobs on pfe."
- "Show the last 30 lines of `/nobackup/<username>/run1/STDOUT.0000` and tell me if the run is still stepping."
- "Copy `plot_fields.py` to `pfe:/nobackup/<username>/analysis/` and run it there."
- "Make a copy of `/nobackup/<username>/run1` called `run1_test` and change `nTimeSteps` to 10 in the copy."
- "Prepare a job script for this run, show it to me, and wait for my OK before submitting it."

Claude turns each request into commands like these and shows you every one before or as it runs it:

```bash
ssh -o BatchMode=yes pfe '/PBS/bin/qstat -u <username>'
ssh -o BatchMode=yes pfe 'tail -30 /nobackup/<username>/run1/STDOUT.0000'
scp -o BatchMode=yes plot_fields.py pfe:/nobackup/<username>/analysis/
```

`BatchMode=yes` means a command fails at once if the connection has closed, instead of waiting for a passcode nobody will type.

### 4. Close the connection when you are done

In Terminal, or after `!` at the Claude prompt:

```bash
ssh -O exit pfe
```

If you forget, it closes by itself 12 hours after the last command that used it.

> **Daily routine in three lines:** `claude` → `! ssh -fN pfe` (passcode + PIN) → ask away. Finish with `ssh -O exit pfe`.

## Waiting for long jobs

Ask "watch job 12345 and tell me when it finishes." Claude starts a small loop on your Mac that checks `qstat` every few minutes and tells you when the job leaves the queue; it then reads the job's output and summarises it. The loop runs on your Mac, not on Pleiades, and each check is one quick command.

Claude Code stops background watchers after about 2 hours. For longer jobs, ask Claude to check again later ("are my jobs done?"), or get the e-mail from PBS by adding `#PBS -m abe` to your job script.

If the connection expires while a job runs, the job is not affected. Run `! ssh -fN pfe` again and carry on.

## Security and data

- Your passcode is typed only by you, into your own terminal, after `!`. Claude cannot log in on its own.
- While the connection is open, anyone who can use your Mac account can reach Pleiades through it, just as they could use your logged-in Terminal. Lock your screen and close the connection at the end of the day.
- Claude reads the output of every command it runs, and that output is sent to Anthropic as part of the conversation. Do not have it print data that your project or NAS policy restricts from leaving the system.
- Submitting jobs uses your allocation. The `CLAUDE.md` rules above make Claude ask first.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Permission denied` or `Host key verification failed` when Claude runs a command | No open connection, so the command refuses to prompt for a passcode | Run `! ssh -fN pfe`, then `ssh -O check pfe` |
| `Control socket connect(...): No such file or directory` | The connection expired or was closed | Run `! ssh -fN pfe` again |
| `! ssh -fN pfe` ends with `Connection closed` | The login did not finish (passcode timed out or mistyped), or an old socket file is in the way | Run `ssh -O exit pfe`, then `rm ~/.ssh/cm-*`, then try again |
| `ssh -O check pfe` says `No ControlPath specified` | The four lines are missing or not at the top of `~/.ssh/config` | Redo setup step 3 |
| `qstat: command not found` | Remote commands run in a short shell without your usual `PATH` | Use full paths (`/PBS/bin/qstat`, `/PBS/bin/qsub`), or wrap the command: `bash -lc '...'` |
| `module: command not found`, or a Python package is missing on `pfe` | Same cause: no login shell, so modules are not loaded | Wrap the command in `bash -lc 'module load ...; ...'` |
| Claude asks for your passcode or tries to log in | It does not know the set-up | Tell it to use the open connection only, and add the `CLAUDE.md` rules above |

## Quick reference

| What | Command |
|---|---|
| Open the shared connection (from Claude Code) | `! ssh -fN pfe` |
| Is it open? | `ssh -O check pfe` |
| Close it | `ssh -O exit pfe` |
| Run one command on Pleiades | `ssh -o BatchMode=yes pfe '<command>'` |
| Copy a file to Pleiades | `scp -o BatchMode=yes file pfe:/nobackup/<username>/` |
| Copy a file back | `scp -o BatchMode=yes pfe:/nobackup/<username>/file .` |
| Your jobs | `ssh -o BatchMode=yes pfe '/PBS/bin/qstat -u <username>'` |
| Undo the setup | `cp ~/.ssh/config.backup ~/.ssh/config` |
