# DrLevy.AI Pilot

> An AI-orchestrated pharmacometrics analysis desktop application — powered by a hierarchy of specialised LLM agents that plan, execute, and document a full population PK/PD modelling workflow from raw data through final regulatory report.

---


DrLevy.AI Pilot converts a pharmacometrician's project description into a fully executed analysis — producing NONMEM control streams, R scripts, diagnostic plots, VPC documentation, and a regulatory-grade final report — all inside a sandboxed Electron desktop window with no cloud storage required.

All session data lives in the browser's **IndexedDB** (via Dexie.js) inside the Electron renderer process. The only network traffic is outbound HTTPS to:

- **OpenRouter API** — LLM inference (free-tier models, no per-user billing)
- **PHP feedback backend** — optional, for collecting user feedback

---

## Key Features

- Six specialised AI agents: Paul (Orchestrator), Project Manager, Data Manager, Pharmacometrician, Medical Writer, QC Manager
- 77-task / 11-section regulatory-grade pharmacometrics checklist
- Human-in-the-loop (HITL) pause-and-review mode
- Fully autonomous execution loop (up to 25 tasks per run)
- Monaco-powered in-app file editor with version history (5 snapshots/file)
- Version diff & one-click revert
- Virtual file system with drag-and-drop dataset upload
- Download All workspace files as ZIP
- Per-user namespaced IndexedDB (multiple analysts can share a machine)
- Hash-based URL routing — each session has a stable URL that survives reload
- Runs 100% offline after first launch (except LLM API calls)

---
